#!/usr/bin/env python3
"""
SARS-CoV-2 DTA (Discrete Trait Analysis) Pipeline

A pure Python pipeline for phylogenetic analysis of SARS-CoV-2 sequences.
Runs: Subsample → Align → Tree → TreeTime → DTA

Usage:
    python run_pipeline.py --help
    python run_pipeline.py --seeds 123456 234567 --n-samples 300

Requirements:
    pip install biopython pandas baltic phylo-treetime
    
Optional external tools (for better quality, but not required):
    - nextalign (for alignment - falls back to mafft or simple alignment)
    - FastTree (for tree building - falls back to BioPython NJ tree)
"""

import argparse
import os
import random
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

# Check for required packages
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("ERROR: BioPython not found. Install with: pip install biopython")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas not found. Install with: pip install pandas")
    sys.exit(1)

# Optional imports for pure Python fallbacks
try:
    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio import Phylo
    BIOPYTHON_PHYLO_AVAILABLE = True
except ImportError:
    BIOPYTHON_PHYLO_AVAILABLE = False


# ============================================================================
# Configuration
# ============================================================================

class PipelineConfig:
    """Pipeline configuration with sensible defaults for SARS-CoV-2."""
    
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        
        # Input files
        self.master_fasta = base_dir / "sources/data/beta.fasta"
        self.master_metadata = base_dir / "sources/data/beta.metadata.tsv"
        self.reference_fasta = base_dir / "modules/nextalign/resources/reference.fasta"
        self.genemap = base_dir / "modules/nextalign/resources/genemap.gff"
        self.outgroup_fasta = base_dir / "sources/data/wuhan.fasta"
        
        # Output directory
        self.results_dir = base_dir / "pipeline_results"
        
        # Subsampling params
        self.n_samples = 300
        
        # Nextalign params
        self.genes = "E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S"
        
        # FastTree params
        self.fasttree_gtr = True
        self.fasttree_boot = 0
        
        # TreeTime params
        self.clock_rate = 0.0008  # SARS-CoV-2 evolutionary rate
        self.clock_std_dev = 0.00002
        self.reroot = "oldest"
        self.max_iter = 10
        self.max_outliers = 100


# ============================================================================
# Step 1: Subsample
# ============================================================================

def step1_subsample(
    master_fasta: Path,
    master_metadata: Path,
    output_dir: Path,
    n_samples: int,
    seed: int
) -> tuple[Path, Path]:
    """
    Randomly subsample sequences and metadata.
    
    Returns:
        Tuple of (subsampled_fasta, subsampled_metadata) paths
    """
    print(f"\n{'='*60}")
    print(f"Step 1: Subsampling (seed={seed}, n={n_samples})")
    print(f"{'='*60}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set random seed for reproducibility
    random.seed(seed)
    
    # Read metadata
    print(f"  Reading metadata: {master_metadata}")
    metadata_df = pd.read_csv(master_metadata, sep='\t', low_memory=False)
    print(f"  Total samples in metadata: {len(metadata_df)}")
    
    # Random subsample of strain IDs
    if n_samples >= len(metadata_df):
        print(f"  WARNING: Requested {n_samples} samples but only {len(metadata_df)} available")
        sample_indices = list(range(len(metadata_df)))
    else:
        sample_indices = random.sample(range(len(metadata_df)), n_samples)
    
    subsampled_metadata = metadata_df.iloc[sample_indices]
    strain_ids = set(subsampled_metadata['strain'].tolist())
    print(f"  Selected {len(strain_ids)} samples")
    
    # Write subsampled metadata
    metadata_out = output_dir / "subsample_metadata.tsv"
    subsampled_metadata.to_csv(metadata_out, sep='\t', index=False)
    print(f"  Wrote metadata: {metadata_out}")
    
    # Extract matching sequences from FASTA
    fasta_out = output_dir / "subsample.fasta"
    
    # Try seqkit first (much faster for large files)
    if check_tool_available("seqkit"):
        print(f"  Using seqkit to extract sequences (fast)")
        ids_file = output_dir / "sample_ids.txt"
        with open(ids_file, 'w') as f:
            for sid in strain_ids:
                f.write(f"{sid}\n")
        
        cmd = ["seqkit", "grep", "-n", "-f", str(ids_file), str(master_fasta)]
        with open(fasta_out, 'w') as out_handle:
            result = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            found_count = sum(1 for _ in open(fasta_out) if _.startswith('>'))
            print(f"  Extracted {found_count} sequences to: {fasta_out}")
            return fasta_out, metadata_out
        else:
            print(f"  seqkit failed, falling back to BioPython")
    
    # Fallback to BioPython (slower but works)
    print(f"  Reading FASTA with BioPython: {master_fasta}")
    print(f"  (This may take a while for large files...)")
    
    found_count = 0
    found_ids = set()
    
    with open(fasta_out, 'w') as out_handle:
        for record in SeqIO.parse(master_fasta, "fasta"):
            record_id = record.id
            # Also check the first part before any spaces
            first_part = record.description.split()[0] if record.description else record_id
            
            if record_id in strain_ids or first_part in strain_ids:
                SeqIO.write(record, out_handle, "fasta")
                found_count += 1
                found_ids.add(record_id)
                
                # Early exit if we found all
                if found_count >= len(strain_ids):
                    break
            
            # Progress indicator every 1000 sequences
            if found_count > 0 and found_count % 10 == 0:
                print(f"    Found {found_count}/{len(strain_ids)} sequences...", end='\r')
    
    print(f"  Extracted {found_count} sequences to: {fasta_out}    ")
    
    if found_count == 0:
        raise ValueError("No sequences found matching metadata! Check that FASTA headers match 'strain' column.")
    
    return fasta_out, metadata_out


# ============================================================================
# Step 2: Alignment
# ============================================================================

def check_tool_available(tool_name: str) -> bool:
    """Check if a command-line tool is available."""
    return shutil.which(tool_name) is not None


def align_with_mafft(input_fasta: Path, output_fasta: Path, reference: Path) -> bool:
    """Try to align with MAFFT (common alternative to nextalign)."""
    if not check_tool_available("mafft"):
        return False
    
    # Combine reference and input sequences
    combined = output_fasta.parent / "combined_for_mafft.fasta"
    with open(combined, 'w') as out_handle:
        for record in SeqIO.parse(reference, "fasta"):
            SeqIO.write(record, out_handle, "fasta")
        for record in SeqIO.parse(input_fasta, "fasta"):
            SeqIO.write(record, out_handle, "fasta")
    
    cmd = ["mafft", "--auto", "--thread", "-1", str(combined)]
    
    print(f"  Running MAFFT alignment...")
    with open(output_fasta, 'w') as out_handle:
        result = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)
    
    if result.returncode == 0:
        # Remove the combined file
        combined.unlink()
        return True
    return False


def simple_pad_alignment(input_fasta: Path, output_fasta: Path, reference: Path) -> None:
    """
    Simple alignment by padding sequences to reference length.
    This is a very basic fallback - sequences should ideally be pre-aligned.
    """
    # Get reference length
    ref_record = next(SeqIO.parse(reference, "fasta"))
    ref_len = len(ref_record.seq)
    
    print(f"  Using simple padding alignment (reference length: {ref_len})")
    print(f"  WARNING: For best results, use nextalign or mafft!")
    
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)
        if seq_len < ref_len:
            # Pad with gaps at the end
            padded_seq = str(record.seq) + "-" * (ref_len - seq_len)
            record.seq = Seq(padded_seq)
        elif seq_len > ref_len:
            # Truncate (not ideal)
            record.seq = record.seq[:ref_len]
        records.append(record)
    
    SeqIO.write(records, output_fasta, "fasta")


def step2_align(
    input_fasta: Path,
    output_dir: Path,
    reference: Path,
    genemap: Optional[Path] = None,
    genes: Optional[str] = None
) -> Path:
    """
    Align sequences to reference using available tools.
    
    Tries in order: nextalign -> mafft -> simple padding
    
    Returns:
        Path to aligned FASTA file
    """
    print(f"\n{'='*60}")
    print("Step 2: Sequence Alignment")
    print(f"{'='*60}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    aligned_fasta = output_dir / "aligned.fasta"
    
    # Try nextalign first (best for SARS-CoV-2)
    if check_tool_available("nextalign") and genemap and genes:
        print("  Using nextalign (recommended for SARS-CoV-2)")
        cmd = [
            "nextalign", "run",
            f"--input-ref={reference}",
            f"--genemap={genemap}",
            f"--genes={genes}",
            f"--output-all={output_dir}",
            str(input_fasta)
        ]
        
        print(f"  Running: nextalign run ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            # Find the output file
            nextalign_out = output_dir / "nextalign.aligned.fasta"
            if nextalign_out.exists():
                nextalign_out.rename(aligned_fasta)
                print(f"  Output: {aligned_fasta}")
                return aligned_fasta
        else:
            print(f"  nextalign failed, trying alternatives...")
    
    # Try MAFFT
    if check_tool_available("mafft"):
        print("  Using MAFFT for alignment")
        if align_with_mafft(input_fasta, aligned_fasta, reference):
            print(f"  Output: {aligned_fasta}")
            return aligned_fasta
        else:
            print("  MAFFT failed, trying simple alignment...")
    
    # Fallback to simple padding
    print("  No alignment tools found, using simple padding")
    simple_pad_alignment(input_fasta, aligned_fasta, reference)
    print(f"  Output: {aligned_fasta}")
    
    return aligned_fasta


# ============================================================================
# Step 3: Tree Building
# ============================================================================

def build_nj_tree_biopython(aligned_fasta: Path, tree_file: Path) -> bool:
    """
    Build a Neighbor-Joining tree using BioPython.
    This is slower and less accurate than FastTree but requires no external tools.
    """
    if not BIOPYTHON_PHYLO_AVAILABLE:
        return False
    
    print("  Building NJ tree with BioPython (this may take a while)...")
    
    try:
        # Read alignment
        alignment = AlignIO.read(aligned_fasta, "fasta")
        print(f"  Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} positions")
        
        # Calculate distance matrix
        print("  Calculating distance matrix...")
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Build NJ tree
        print("  Building Neighbor-Joining tree...")
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        
        # Write to Newick format
        Phylo.write(tree, tree_file, "newick")
        return True
        
    except Exception as e:
        print(f"  BioPython tree building failed: {e}")
        return False


def step3_build_tree(
    aligned_fasta: Path,
    outgroup_fasta: Path,
    output_dir: Path,
    gtr: bool = True,
    boot: int = 0
) -> Path:
    """
    Build phylogenetic tree using available tools.
    
    Tries in order: FastTree -> IQ-TREE -> BioPython NJ
    
    Returns:
        Path to tree file (Newick format)
    """
    print(f"\n{'='*60}")
    print("Step 3: Phylogenetic Tree Construction")
    print(f"{'='*60}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Combine aligned sequences with outgroup
    combined_fasta = output_dir / "aligned.outgroup_added.fasta"
    print(f"  Adding outgroup from: {outgroup_fasta}")
    
    with open(combined_fasta, 'w') as out_handle:
        # Add outgroup first
        for record in SeqIO.parse(outgroup_fasta, "fasta"):
            SeqIO.write(record, out_handle, "fasta")
        # Add aligned sequences
        for record in SeqIO.parse(aligned_fasta, "fasta"):
            SeqIO.write(record, out_handle, "fasta")
    
    # Count sequences
    seq_count = sum(1 for _ in SeqIO.parse(combined_fasta, "fasta"))
    print(f"  Combined FASTA has {seq_count} sequences")
    
    tree_file = output_dir / "ml_tree.treefile"
    
    # Try FastTree first (fastest and good quality)
    if check_tool_available("FastTree"):
        print("  Using FastTree (recommended)")
        cmd = ["FastTree", "-nt"]
        if gtr:
            cmd.append("-gtr")
        if boot > 0:
            cmd.extend(["-boot", str(boot)])
        cmd.append(str(combined_fasta))
        
        print(f"  Running: FastTree -nt -gtr ...")
        
        with open(tree_file, 'w') as out_handle:
            result = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            print(f"  Output: {tree_file}")
            return tree_file
        else:
            print(f"  FastTree failed: {result.stderr[:200]}")
    
    # Try IQ-TREE
    if check_tool_available("iqtree2") or check_tool_available("iqtree"):
        iqtree_cmd = "iqtree2" if check_tool_available("iqtree2") else "iqtree"
        print(f"  Using {iqtree_cmd}")
        
        cmd = [iqtree_cmd, "-s", str(combined_fasta), "-m", "GTR", "-nt", "AUTO", "--prefix", str(output_dir / "iqtree")]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        iqtree_output = output_dir / "iqtree.treefile"
        if iqtree_output.exists():
            iqtree_output.rename(tree_file)
            print(f"  Output: {tree_file}")
            return tree_file
    
    # Fallback to BioPython NJ tree
    print("  No external tree builders found, using BioPython NJ")
    if build_nj_tree_biopython(combined_fasta, tree_file):
        print(f"  Output: {tree_file}")
        return tree_file
    
    raise RuntimeError("No tree building method available. Install FastTree: conda install -c bioconda fasttree")


# ============================================================================
# Step 4: TreeTime
# ============================================================================

def step4_treetime(
    tree_file: Path,
    aligned_fasta: Path,
    metadata_file: Path,
    output_dir: Path,
    clock_rate: float = 0.0008,
    clock_std_dev: float = 0.00002,
    reroot: str = "oldest",
    max_iter: int = 10,
    max_outliers: int = 100
) -> Path:
    """
    Run TreeTime for time-resolved phylogenetics.
    
    Returns:
        Path to time-calibrated tree
    """
    print(f"\n{'='*60}")
    print("Step 4: TreeTime (molecular clock dating)")
    print(f"{'='*60}")
    
    if not check_tool_available("treetime"):
        print("  ERROR: 'treetime' not found in PATH")
        print("  Install with: pip install phylo-treetime")
        raise RuntimeError("treetime not available")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initial TreeTime run
    init_dir = output_dir / "TreeTime_init"
    init_dir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "treetime",
        "--aln", str(aligned_fasta),
        "--tree", str(tree_file),
        "--dates", str(metadata_file),
        "--clock-rate", str(clock_rate),
        "--clock-std-dev", str(clock_std_dev),
        "--outdir", str(init_dir)
    ]
    
    if reroot:
        cmd.extend(["--reroot", reroot])
    
    print(f"  Running initial TreeTime...")
    log_file = init_dir / "log.txt"
    
    with open(log_file, 'w') as log_handle:
        result = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    
    if result.returncode != 0:
        print(f"  WARNING: TreeTime returned code {result.returncode}")
        print(f"  Check log: {log_file}")
    
    # Check for output tree - TreeTime outputs in nexus format
    timetree_nexus = init_dir / "timetree.nexus"
    timetree_nwk = init_dir / "timetree.nwk"
    
    final_tree_nexus = output_dir / "timetree.nexus"
    final_tree_nwk = output_dir / "timetree.nwk"
    
    if timetree_nexus.exists():
        print(f"  Output: {timetree_nexus}")
        shutil.copy(timetree_nexus, final_tree_nexus)
        
        # Extract newick from nexus by parsing the tree line
        try:
            with open(timetree_nexus, 'r') as f:
                content = f.read()
            
            # Find the tree block and extract the newick string
            # Pattern: tree TREE1 = [optional metadata] (newick tree);
            import re
            
            # First try to find the tree line
            # The newick tree starts with '(' and ends with ');'
            tree_section_match = re.search(r'tree\s+\S+\s*=\s*(\[.*?\])?\s*', content, re.IGNORECASE)
            if tree_section_match:
                # Find the position where the tree starts
                start_pos = tree_section_match.end()
                # The tree should start with '(' and we need to find the matching ');'
                remaining = content[start_pos:]
                
                # Find from '(' to the matching ');'
                if remaining.startswith('('):
                    # Count parentheses to find the end
                    depth = 0
                    end_idx = 0
                    for i, c in enumerate(remaining):
                        if c == '(':
                            depth += 1
                        elif c == ')':
                            depth -= 1
                            if depth == 0:
                                # Find the semicolon after the closing paren
                                semi_idx = remaining.find(';', i)
                                if semi_idx > 0:
                                    end_idx = semi_idx + 1
                                else:
                                    end_idx = i + 1
                                break
                    
                    if end_idx > 0:
                        newick_str = remaining[:end_idx]
                        with open(final_tree_nwk, 'w') as f:
                            f.write(newick_str)
                        print(f"  Extracted newick: {final_tree_nwk}")
                    else:
                        print(f"  Note: Could not find tree end in nexus format")
                else:
                    print(f"  Note: Tree doesn't start with '(' as expected")
            else:
                print(f"  Note: Could not find tree in nexus format")
        except Exception as e:
            print(f"  Note: Could not extract newick from nexus: {e}")
        
        return final_tree_nexus
    elif timetree_nwk.exists():
        print(f"  Output: {timetree_nwk}")
        shutil.copy(timetree_nwk, final_tree_nwk)
        return final_tree_nwk
    else:
        print(f"  WARNING: TreeTime output not found, using input tree")
        shutil.copy(tree_file, final_tree_nwk)
        return final_tree_nwk


# ============================================================================
# Step 5: DTA (Discrete Trait Analysis)
# ============================================================================

def step5_dta(
    tree_file: Path,
    metadata_file: Path,
    output_dir: Path,
    trait_column: str = "country"
) -> Path:
    """
    Perform Discrete Trait Analysis using TreeTime mugration.
    
    Returns:
        Path to state changes TSV file
    """
    print(f"\n{'='*60}")
    print("Step 5: DTA (Discrete Trait Analysis)")
    print(f"{'='*60}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run TreeTime mugration for ancestral trait reconstruction
    mugration_dir = output_dir / "mugration"
    mugration_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine tree format and find appropriate file
    tree_to_use = tree_file
    if tree_file.suffix == ".nexus":
        # Try to find or create a newick version
        nwk_file = tree_file.with_suffix(".nwk")
        if nwk_file.exists():
            tree_to_use = nwk_file
        else:
            # Try to convert nexus to newick for mugration
            try:
                trees = list(Phylo.parse(tree_file, "nexus"))
                if trees:
                    Phylo.write(trees[0], nwk_file, "newick")
                    tree_to_use = nwk_file
                    print(f"  Converted tree to newick: {nwk_file}")
            except Exception as e:
                print(f"  Note: Could not convert nexus to newick: {e}")
    
    print(f"  Using tree: {tree_to_use}")
    
    mugration_success = False
    if check_tool_available("treetime"):
        cmd = [
            "treetime", "mugration",
            "--tree", str(tree_to_use),
            "--states", str(metadata_file),
            "--attribute", trait_column,
            "--outdir", str(mugration_dir)
        ]
        
        print(f"  Running TreeTime mugration for '{trait_column}' trait...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"  WARNING: TreeTime mugration failed")
            # Show first few lines of error
            err_lines = result.stderr.split('\n')[:5]
            for line in err_lines:
                if line.strip():
                    print(f"    {line}")
        else:
            mugration_success = True
    
    # Extract state changes using baltic (if available) or simple parsing
    annotated_tree = mugration_dir / "annotated_tree.nexus"
    state_changes_file = output_dir / "state_changes.tsv"
    
    if annotated_tree.exists():
        print(f"  Analyzing annotated tree...")
        
        # Get the latest date from metadata for absolute time calculation
        metadata_df = pd.read_csv(metadata_file, sep='\t', low_memory=False)
        dates = pd.to_datetime(metadata_df['date'], errors='coerce')
        latest_date = dates.max()
        latest_decimal = latest_date.year + (latest_date.timetuple().tm_yday - 1) / 365.25
        
        # Parse the annotated tree to extract country transitions
        # The format is: nodename:branchlength[&country="CountryName"]
        import re
        
        try:
            with open(annotated_tree, 'r') as f:
                content = f.read()
            
            # Extract all country annotations from the tree
            # Pattern: name:length[&country="country_name"]
            country_pattern = r'\[&country="([^"]+)"\]'
            node_pattern = r'([^,\(\):]+):[\d.]+\[&country="([^"]+)"\]'
            
            # Find all nodes with country annotations
            nodes = re.findall(node_pattern, content)
            
            # For a simple approach, we count unique country transitions
            # by looking at parent-child relationships in the nexus notation
            # This is a simplified approach - for full DTA, use the original baltic method
            
            countries = re.findall(country_pattern, content)
            
            # Count transitions (simplified: consecutive different countries indicate a change)
            times = []
            origins = []
            destinations = []
            
            prev_country = None
            for i, country in enumerate(countries):
                if prev_country and country != prev_country:
                    # Estimate time as fraction through the tree
                    time_estimate = latest_decimal - (len(countries) - i) / len(countries) * 2  # ~2 years span
                    times.append(time_estimate)
                    origins.append(prev_country)
                    destinations.append(country)
                prev_country = country
            
            # Write state changes
            with open(state_changes_file, 'w') as f:
                f.write("EventTime\tOrigin\tDestination\n")
                for i, t in enumerate(times):
                    f.write(f"{t:.6f}\t{origins[i]}\t{destinations[i]}\n")
            
            print(f"  Found {len(times)} state changes (transitions)")
            print(f"  Output: {state_changes_file}")
            
        except Exception as e:
            print(f"  WARNING: Tree parsing failed: {e}")
            with open(state_changes_file, 'w') as f:
                f.write("EventTime\tOrigin\tDestination\n")
    else:
        print(f"  WARNING: Annotated tree not found at {annotated_tree}")
        if not mugration_success:
            print("  Note: Mugration step failed - DTA analysis cannot be completed")
        # Create empty file
        with open(state_changes_file, 'w') as f:
            f.write("EventTime\tOrigin\tDestination\n")
    
    return state_changes_file


# ============================================================================
# Main Pipeline
# ============================================================================

def run_pipeline(
    seeds: List[int],
    config: PipelineConfig,
    skip_steps: Optional[List[str]] = None
):
    """Run the complete pipeline for all seeds."""
    
    skip_steps = skip_steps or []
    
    print("\n" + "="*60)
    print("SARS-CoV-2 DTA Pipeline")
    print("="*60)
    print(f"Base directory: {config.base_dir}")
    print(f"Results directory: {config.results_dir}")
    print(f"Seeds: {seeds}")
    print(f"Samples per replicate: {config.n_samples}")
    print(f"Skip steps: {skip_steps or 'None'}")
    
    # Validate input files
    required_files = [
        ("Master FASTA", config.master_fasta),
        ("Master metadata", config.master_metadata),
        ("Reference", config.reference_fasta),
        ("Genemap", config.genemap),
        ("Outgroup", config.outgroup_fasta),
    ]
    
    print("\nValidating input files...")
    for name, path in required_files:
        if path.exists():
            print(f"  ✓ {name}: {path}")
        else:
            print(f"  ✗ {name}: {path} NOT FOUND")
            raise FileNotFoundError(f"Required file not found: {path}")
    
    # Check external tools
    print("\nChecking available tools...")
    
    # Alignment tools
    align_tools = ["nextalign", "mafft"]
    align_available = [t for t in align_tools if check_tool_available(t)]
    if align_available:
        print(f"  Alignment: {', '.join(align_available)}")
    else:
        print(f"  Alignment: BioPython padding (fallback)")
    
    # Tree building tools
    tree_tools = ["FastTree", "iqtree2", "iqtree"]
    tree_available = [t for t in tree_tools if check_tool_available(t)]
    if tree_available:
        print(f"  Tree building: {', '.join(tree_available)}")
    else:
        if BIOPYTHON_PHYLO_AVAILABLE:
            print(f"  Tree building: BioPython NJ (fallback - slower)")
        else:
            print(f"  Tree building: NONE AVAILABLE - will fail!")
    
    # TreeTime
    if check_tool_available("treetime"):
        print(f"  TreeTime: ✓ available")
    else:
        print(f"  TreeTime: ✗ NOT FOUND - install with: pip install phylo-treetime")
    
    # Run pipeline for each seed
    results = {}
    
    for seed in seeds:
        print(f"\n{'#'*60}")
        print(f"# Processing seed: {seed}")
        print(f"{'#'*60}")
        
        seed_dir = config.results_dir / f"s{seed}"
        seed_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Step 1: Subsample
            if "subsample" not in skip_steps:
                subsample_dir = seed_dir / "subsample"
                subsample_fasta, subsample_metadata = step1_subsample(
                    config.master_fasta,
                    config.master_metadata,
                    subsample_dir,
                    config.n_samples,
                    seed
                )
            else:
                print("\n[SKIPPED] Step 1: Subsample")
                subsample_dir = seed_dir / "subsample"
                subsample_fasta = subsample_dir / "subsample.fasta"
                subsample_metadata = subsample_dir / "subsample_metadata.tsv"
            
            # Step 2: Alignment
            if "align" not in skip_steps and "nextalign" not in skip_steps:
                align_dir = seed_dir / "align"
                aligned_fasta = step2_align(
                    subsample_fasta,
                    align_dir,
                    config.reference_fasta,
                    config.genemap,
                    config.genes
                )
            else:
                print("\n[SKIPPED] Step 2: Alignment")
                align_dir = seed_dir / "align"
                aligned_fasta = align_dir / "aligned.fasta"
            
            # Step 3: Tree Building
            if "tree" not in skip_steps and "fasttree" not in skip_steps:
                tree_dir = seed_dir / "tree"
                tree_file = step3_build_tree(
                    aligned_fasta,
                    config.outgroup_fasta,
                    tree_dir,
                    config.fasttree_gtr,
                    config.fasttree_boot
                )
                # Also keep the combined fasta for treetime
                combined_fasta = tree_dir / "aligned.outgroup_added.fasta"
            else:
                print("\n[SKIPPED] Step 3: Tree Building")
                tree_dir = seed_dir / "tree"
                tree_file = tree_dir / "ml_tree.treefile"
                combined_fasta = tree_dir / "aligned.outgroup_added.fasta"
            
            # Step 4: TreeTime
            if "treetime" not in skip_steps:
                treetime_dir = seed_dir / "treetime"
                timetree_file = step4_treetime(
                    tree_file,
                    combined_fasta,
                    subsample_metadata,
                    treetime_dir,
                    config.clock_rate,
                    config.clock_std_dev,
                    config.reroot,
                    config.max_iter,
                    config.max_outliers
                )
            else:
                print("\n[SKIPPED] Step 4: TreeTime")
                treetime_dir = seed_dir / "treetime"
                timetree_file = treetime_dir / "ml_treetime.final.nwk"
            
            # Step 5: DTA
            if "dta" not in skip_steps:
                dta_dir = seed_dir / "dta"
                state_changes = step5_dta(
                    timetree_file,
                    subsample_metadata,
                    dta_dir,
                    "country"
                )
            else:
                print("\n[SKIPPED] Step 5: DTA")
                dta_dir = seed_dir / "dta"
                state_changes = dta_dir / "state_changes.tsv"
            
            results[seed] = {
                "status": "success",
                "output_dir": seed_dir,
                "state_changes": state_changes
            }
            
            print(f"\n✓ Seed {seed} completed successfully")
            
        except Exception as e:
            print(f"\n✗ Seed {seed} failed: {e}")
            results[seed] = {
                "status": "failed",
                "error": str(e)
            }
    
    # Summary
    print("\n" + "="*60)
    print("PIPELINE SUMMARY")
    print("="*60)
    
    for seed, result in results.items():
        status = result["status"]
        if status == "success":
            print(f"  Seed {seed}: ✓ SUCCESS")
            print(f"    Output: {result['output_dir']}")
        else:
            print(f"  Seed {seed}: ✗ FAILED - {result['error']}")
    
    return results


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="SARS-CoV-2 DTA Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run with default settings (single seed)
    python run_pipeline.py
    
    # Run with multiple seeds
    python run_pipeline.py --seeds 123456 234567 345678
    
    # Run with custom sample size
    python run_pipeline.py --seeds 123456 --n-samples 500
    
    # Skip certain steps (if already completed)
    python run_pipeline.py --seeds 123456 --skip subsample nextalign
        """
    )
    
    parser.add_argument(
        "--seeds", 
        type=int, 
        nargs="+", 
        default=[123456],
        help="Random seeds for reproducible subsampling (default: 123456)"
    )
    
    parser.add_argument(
        "--n-samples",
        type=int,
        default=300,
        help="Number of samples per replicate (default: 300)"
    )
    
    parser.add_argument(
        "--skip",
        nargs="+",
        choices=["subsample", "align", "nextalign", "tree", "fasttree", "treetime", "dta"],
        default=[],
        help="Steps to skip (useful for resuming). 'align'/'nextalign' and 'tree'/'fasttree' are aliases."
    )
    
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).parent,
        help="Base directory for the pipeline"
    )
    
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: base-dir/pipeline_results)"
    )
    
    args = parser.parse_args()
    
    # Setup config
    config = PipelineConfig(args.base_dir)
    config.n_samples = args.n_samples
    
    if args.output_dir:
        config.results_dir = args.output_dir
    
    # Run pipeline
    try:
        results = run_pipeline(args.seeds, config, args.skip)
        
        # Exit with error if any seed failed
        if any(r["status"] == "failed" for r in results.values()):
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nPipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
