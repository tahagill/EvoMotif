# EvoMotif: Complete Technical Guide

**Version**: 0.1.0  
**Last Updated**: December 23, 2025

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Algorithm Details](#algorithm-details)
4. [Statistical Methods](#statistical-methods)
5. [API Reference](#api-reference)
6. [Implementation Guide](#implementation-guide)
7. [Use Cases](#use-cases)
8. [Limitations](#limitations)
9. [Troubleshooting](#troubleshooting)

---

## Overview

EvoMotif discovers evolutionarily conserved protein motifs by analyzing multiple sequence alignments. The pipeline combines information theory (Shannon entropy), evolutionary substitution matrices (BLOSUM62), permutation-based statistical testing, and FDR correction to identify functionally important regions with high confidence.

**Core principle**: If a residue or motif is conserved across millions of years of evolution, it likely serves an important structural or functional role.

### Pipeline Stages

1. **Retrieval**: Query NCBI for homologous sequences
2. **Alignment**: Multiple sequence alignment with MAFFT
3. **Conservation**: Calculate per-position conservation scores
4. **Discovery**: Identify conserved motifs using sliding windows
5. **Validation**: Statistical testing with permutation and FDR
6. **Phylogeny**: Build evolutionary tree with FastTree
7. **Structure**: Map conservation to 3D coordinates (optional)

---

## Installation

### System Requirements

- Python 3.8 or higher
- 8GB RAM minimum (16GB recommended for large proteins)
- Linux, macOS, or Windows with WSL
- Internet connection (for NCBI queries)

### External Dependencies

EvoMotif requires three bioinformatics tools:

**MAFFT** - Multiple sequence alignment
```bash
# Ubuntu/Debian
sudo apt-get install mafft

# macOS
brew install mafft

# Conda
conda install -c bioconda mafft

# Verify installation
mafft --version  # Should show v7.x or higher
```

**FastTree** - Phylogenetic tree inference
```bash
# Ubuntu/Debian
sudo apt-get install fasttree

# macOS
brew install fasttree

# Conda
conda install -c bioconda fasttree

# Verify
FastTree -help  # Should show usage info
```

**DSSP** (Optional) - Secondary structure assignment
```bash
# Ubuntu/Debian
sudo apt-get install dssp

# macOS
brew install dssp

# Conda
conda install -c salilab dssp
```

### Python Installation

```bash
# Clone repository
git clone https://github.com/yourusername/evomotif.git
cd evomotif

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Install
pip install -e .

# Verify
python -c "import evomotif; print('Success')"
```

### Verify Dependencies

```python
import evomotif

pipeline = evomotif.EvoMotifPipeline()
deps = pipeline.check_dependencies()

for tool, available in deps.items():
    status = "✓" if available else "✗"
    print(f"{status} {tool}")
```

Expected output:
```
✓ python
✓ mafft
✓ fasttree
✓ biopython
✓ numpy
✓ scipy
```

---

## Algorithm Details

### 1. Sequence Retrieval

**Purpose**: Obtain diverse, high-quality homologous sequences from NCBI

**Implementation**:

```python
from Bio import Entrez
from evomotif.retrieval import SequenceRetriever

retriever = SequenceRetriever(email="your@email.com")
sequences = retriever.fetch_sequences("hemoglobin", max_results=100)
```

**Algorithm steps**:

1. **Query NCBI Protein database** with protein name
2. **Fetch sequence records** in batches (100 per batch for efficiency)
3. **Quality filtering**:
   - Remove sequences with ambiguous residues (X, B, Z > 5% of length)
   - Remove sequences shorter than 50% of median length
   - Remove sequences with >80% low complexity regions
4. **Redundancy reduction**:
   - Cluster sequences by 95% identity threshold
   - Select representative from each cluster (longest sequence)
5. **Taxonomic diversity**:
   - Prioritize sequences from different orders/classes
   - Prefer model organisms with well-annotated genomes

**Why this approach?**
- Diversity increases signal-to-noise ratio for conservation
- Quality control prevents alignment artifacts
- Redundancy reduction speeds computation without losing information

**Parameters**:
- `max_results`: More sequences improve statistics but increase runtime
  - Minimum: 10 (basic analysis)
  - Recommended: 50-100 (robust statistics)
  - Maximum: 500 (comprehensive but slow)

**Output**: List of BioPython SeqRecord objects

### 2. Multiple Sequence Alignment

**Purpose**: Align homologous positions across sequences

**Implementation**:

```python
from evomotif.alignment import SequenceAligner

aligner = SequenceAligner(threads=8)
alignment = aligner.align(sequences, output="aligned.fasta")
```

**Tool**: MAFFT (Multiple Alignment using Fast Fourier Transform)

**Algorithm**: MAFFT uses FFT to detect homologous regions, then applies iterative refinement

**Parameters**:
- `--auto`: Automatic strategy selection based on input size
- `--thread N`: Parallel processing
- `--reorder`: Order sequences by similarity (improves visualization)

**Quality metrics**:

After alignment, check:
```python
from evomotif.alignment import AlignmentQuality

quality = AlignmentQuality(alignment)
print(f"Gap percentage: {quality.gap_percentage():.1f}%")  # <20% is good
print(f"Conservation: {quality.mean_conservation():.3f}")   # >0.5 is good
print(f"Identity: {quality.pairwise_identity():.3f}")       # 40-60% typical
```

**Interpretation**:
- High gaps (>30%): May indicate distant homologs or alignment errors
- Low conservation (<0.4): Protein family may be too divergent
- Very high identity (>90%): Redundancy, consider reducing sequences

### 3. Conservation Scoring

**Purpose**: Quantify evolutionary constraint at each alignment position

EvoMotif combines two complementary metrics:

#### A. Shannon Entropy

**Mathematical definition**:

```
H(i) = -Σ p_a(i) × log₂(p_a(i))
```

where:
- `i` = alignment column (position)
- `a` = amino acid type (20 possibilities)
- `p_a(i)` = frequency of amino acid `a` at position `i`
- `H(i)` = entropy in bits (0 to log₂(20) = 4.32 bits)

**Normalized conservation score**:

```
C_shannon(i) = 1 - H(i) / log₂(20)
```

Range: [0, 1] where 1 = perfect conservation, 0 = maximum variability

**Example calculation**:

```
Position 45: [H, H, H, H, H, H, H, H, H, R]
             (9 His, 1 Arg out of 10 sequences)

p_H = 9/10 = 0.9
p_R = 1/10 = 0.1

H(45) = -(0.9 × log₂(0.9) + 0.1 × log₂(0.1))
      = -(0.9 × -0.152 + 0.1 × -3.322)
      = 0.469 bits

C_shannon(45) = 1 - 0.469/4.322 = 0.891
```

**Interpretation**: Position 45 is highly conserved (89.1%)

**Why Shannon entropy?**
- Grounded in information theory
- Quantifies "surprise" or uncertainty
- High entropy = high uncertainty = low conservation
- No biological assumptions required
- Sensitive to rare variants

#### B. BLOSUM62 Score

**Mathematical definition**:

```
B(i) = (1 / (n(n-1))) × Σ_{j<k} BLOSUM62(a_j, a_k)
```

where:
- `n` = number of sequences
- `a_j, a_k` = amino acids at position `i` in sequences `j` and `k`
- BLOSUM62 = substitution score from empirical protein evolution data

**BLOSUM62 matrix** (excerpt):

```
     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A    4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R   -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
H   -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
...
```

Positive scores = functionally similar  
Negative scores = functionally dissimilar

**Normalized score**:

```
B_norm(i) = (B(i) - B_min) / (B_max - B_min)
```

Scaled to [0, 1]

**Example**:

```
Position 45: [H, H, H, H, H, H, H, H, H, R]

Pairwise comparisons:
- BLOSUM62(H, H) = 8 (identical, 9 × 8 = 72 comparisons)
- BLOSUM62(H, R) = 0 (conservative, 9 comparisons)
- BLOSUM62(R, R) = 5 (would be if more Rs present)

Total pairwise score = 72 × 8 + 9 × 0 = 576
Number of pairs = 10 × 9 / 2 = 45

B(45) = 576 / 45 = 12.8

Normalized (assuming B_min=-4, B_max=15):
B_norm(45) = (12.8 - (-4)) / (15 - (-4)) = 0.884
```

**Why BLOSUM62?**
- Based on real evolutionary data (aligned protein blocks)
- Captures functional constraints beyond identity
- Distinguishes conservative (Leu↔Ile) vs. radical (Lys↔Asp) substitutions
- Widely validated and used in bioinformatics

#### C. Combined Conservation Score

**Formula**:

```
C_final(i) = α × C_shannon(i) + (1-α) × B_norm(i)
```

Default: α = 0.5 (equal weighting)

**Rationale**:
- Shannon: Detects strict conservation (catalytic residues, structural cores)
- BLOSUM: Detects functional conservation (physicochemical similarity)
- Combined: Comprehensive assessment of evolutionary constraint

**Implementation**:

```python
from evomotif.conservation import ConservationScorer

scorer = ConservationScorer()
conservation = scorer.calculate_conservation_scores(
    alignment,
    method="combined",
    weights=(0.5, 0.5)  # (Shannon, BLOSUM)
)
```

**Tuning weights**:
- Higher Shannon weight (0.7, 0.3): Emphasize strict identity
- Higher BLOSUM weight (0.3, 0.7): Emphasize functional equivalence
- Application-specific: Structural studies may prefer Shannon, functional studies BLOSUM

### 4. Motif Discovery

**Purpose**: Identify contiguous regions of high conservation

**Algorithm**: Sliding window with adaptive thresholds

#### Step-by-step process:

**Step 1: Window scanning**

```python
# Pseudocode
for window_size in [7, 9, 11, 13, 15]:
    for start_pos in range(alignment_length - window_size + 1):
        end_pos = start_pos + window_size
        window_conservation = conservation[start_pos:end_pos]
        
        # Calculate window statistics
        mean_cons = np.mean(window_conservation)
        std_cons = np.std(window_conservation)
        gap_freq = calculate_gap_frequency(alignment[:, start_pos:end_pos])
        
        # Apply filters
        if mean_cons >= min_conservation and \
           std_cons <= max_std and \
           gap_freq <= max_gap:
            candidates.append({
                'start': start_pos,
                'end': end_pos,
                'conservation': mean_cons,
                'size': window_size
            })
```

**Step 2: Overlap resolution**

When multiple windows overlap, keep the highest-scoring:

```python
def resolve_overlaps(candidates):
    # Sort by conservation score (descending)
    candidates.sort(key=lambda x: x['conservation'], reverse=True)
    
    selected = []
    for candidate in candidates:
        # Check if overlaps with already selected motifs
        overlaps = False
        for selected_motif in selected:
            if ranges_overlap(candidate, selected_motif):
                overlaps = True
                break
        
        if not overlaps:
            selected.append(candidate)
    
    return selected
```

**Step 3: Consecutive merging**

Merge adjacent high-conservation regions:

```python
def merge_consecutive(motifs):
    merged = []
    current = motifs[0]
    
    for next_motif in motifs[1:]:
        # If motifs are adjacent or overlap slightly
        if next_motif['start'] - current['end'] <= 3:
            # Merge
            current['end'] = next_motif['end']
            current['conservation'] = np.mean([
                current['conservation'],
                next_motif['conservation']
            ])
        else:
            merged.append(current)
            current = next_motif
    
    merged.append(current)
    return merged
```

**Parameters**:

```python
# Default configuration
motif_config = {
    'window_sizes': [7, 9, 11, 13, 15],  # Try multiple sizes
    'min_conservation': 0.70,            # Threshold (adjustable)
    'max_std': 0.15,                     # Consistency within window
    'max_gap': 0.30,                     # Require 70% coverage
    'merge_distance': 3                  # Max gap for merging
}
```

**Tuning guide**:

| Parameter | Lower value | Higher value |
|-----------|-------------|--------------|
| `min_conservation` | More motifs, lower confidence | Fewer motifs, higher confidence |
| `max_std` | Allow variable windows | Require uniform conservation |
| `max_gap` | Allow gappy regions | Require well-aligned regions |
| `window_sizes` | Shorter motifs | Longer domains |

**Example**:

```python
from evomotif.motif_discovery import MotifFinder

finder = MotifFinder()
motifs = finder.find_motifs(
    conservation_scores,
    threshold=0.75,  # Strict
    window_sizes=[7, 9, 11]  # Shorter motifs
)

# Results
for motif in motifs:
    print(f"Motif {motif['sequence']} at {motif['start']}-{motif['end']}")
    print(f"  Conservation: {motif['conservation']:.3f}")
```

---

## Statistical Methods

### 1. Permutation Testing

**Null hypothesis**: Observed conservation pattern is indistinguishable from random

**Test procedure**:

```python
def permutation_test(observed_motif, all_conservation, n_perm=10000):
    """
    Calculate p-value for motif conservation.
    
    Args:
        observed_motif: Dict with 'start', 'end', 'conservation'
        all_conservation: Array of conservation scores for entire alignment
        n_perm: Number of random permutations
    
    Returns:
        p_value: Fraction of random samples >= observed
    """
    # Extract observed statistic
    motif_length = observed_motif['end'] - observed_motif['start']
    observed_score = observed_motif['conservation']
    
    # Generate null distribution
    null_scores = []
    np.random.seed(42)  # Reproducibility
    
    for _ in range(n_perm):
        # Shuffle conservation scores
        shuffled = np.random.permutation(all_conservation)
        
        # Extract random window of same length
        random_start = np.random.randint(0, len(shuffled) - motif_length)
        random_window = shuffled[random_start:random_start + motif_length]
        random_score = np.mean(random_window)
        
        null_scores.append(random_score)
    
    # Calculate p-value
    p_value = np.mean(np.array(null_scores) >= observed_score)
    
    return p_value
```

**Interpretation**:
- p < 0.001: Very strong evidence for conservation
- p < 0.01: Strong evidence
- p < 0.05: Moderate evidence
- p ≥ 0.05: Insufficient evidence

**Why permutation over parametric tests?**
- No distributional assumptions
- Exact p-values for given data
- Appropriate for complex dependencies in sequence data

### 2. Multiple Testing Correction

**Problem**: Testing many motifs increases false positive rate

Example: Testing 50 motifs at α = 0.05, expect ~2.5 false positives by chance

**Solution**: False Discovery Rate (FDR) control via Benjamini-Hochberg

**Algorithm**:

```python
def benjamini_hochberg(p_values, alpha=0.05):
    """
    Control false discovery rate.
    
    Args:
        p_values: Array of p-values from multiple tests
        alpha: Desired FDR level
    
    Returns:
        reject: Boolean array indicating which tests reject H0
        adjusted_p: FDR-adjusted p-values
    """
    m = len(p_values)
    
    # Sort p-values and track original indices
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    # Find largest i where p_i <= (i/m) × alpha
    thresholds = np.arange(1, m + 1) / m * alpha
    reject = sorted_p <= thresholds
    
    if reject.any():
        max_i = np.where(reject)[0][-1]
        # Reject all hypotheses up to max_i
        reject[:max_i + 1] = True
    
    # Calculate adjusted p-values (q-values)
    adjusted_p = np.minimum.accumulate(sorted_p[::-1] * m / np.arange(m, 0, -1))[::-1]
    adjusted_p = np.minimum(adjusted_p, 1.0)
    
    # Restore original order
    original_order = np.argsort(sorted_indices)
    reject = reject[original_order]
    adjusted_p = adjusted_p[original_order]
    
    return reject, adjusted_p
```

**Example**:

```python
# Original p-values for 5 motifs
p_values = [0.001, 0.012, 0.045, 0.078, 0.150]

# Apply BH correction
reject, q_values = benjamini_hochberg(p_values, alpha=0.05)

print("Motif  P-value  Q-value  Reject?")
for i, (p, q, r) in enumerate(zip(p_values, q_values, reject), 1):
    print(f"  {i}    {p:.3f}   {q:.3f}   {'Yes' if r else 'No'}")
```

Output:
```
Motif  P-value  Q-value  Reject?
  1    0.001   0.005   Yes
  2    0.012   0.030   Yes
  3    0.045   0.075   No
  4    0.078   0.098   No
  5    0.150   0.150   No
```

**Why FDR instead of Bonferroni?**
- Bonferroni: α_adj = α / m (very conservative, many false negatives)
- FDR: Controls expected proportion of false discoveries (more powerful)
- Appropriate when multiple true positives expected (common in motif discovery)

### 3. Effect Size

**Purpose**: Quantify practical significance beyond statistical significance

Small p-values can result from large sample sizes even with tiny effects. Effect size measures magnitude.

**Cohen's d**:

```
d = (μ_motif - μ_background) / σ_pooled
```

where:
- μ_motif = mean conservation in motif
- μ_background = mean conservation in non-motif regions
- σ_pooled = pooled standard deviation

**Interpretation**:
- d < 0.2: Small effect (may not be biologically meaningful)
- d = 0.5: Medium effect (noticeable)
- d ≥ 0.8: Large effect (substantial conservation)

**Implementation**:

```python
def cohens_d(motif_scores, background_scores):
    """Calculate Cohen's d effect size."""
    n1 = len(motif_scores)
    n2 = len(background_scores)
    
    mean1 = np.mean(motif_scores)
    mean2 = np.mean(background_scores)
    
    var1 = np.var(motif_scores, ddof=1)
    var2 = np.var(background_scores, ddof=1)
    
    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    d = (mean1 - mean2) / pooled_std
    return d
```

**EvoMotif requirement**: Report motifs only if p < 0.05 (FDR-corrected) AND d > 0.5

---

## API Reference

### High-Level Function

#### `evomotif.analyze_protein()`

Complete protein analysis pipeline in one function.

**Signature**:
```python
evomotif.analyze_protein(
    protein_name: str,
    email: str,
    output_dir: str = None,
    pdb_id: str = None,
    max_sequences: int = 50,
    min_conservation: float = 0.70,
    window_sizes: List[int] = [7, 9, 11, 13, 15],
    threads: int = 4,
    verbose: bool = False
) -> AnalysisResults
```

**Parameters**:
- `protein_name` (str): Protein name for NCBI query (e.g., "hemoglobin", "p53")
- `email` (str): Email for NCBI Entrez (required by NCBI policy)
- `output_dir` (str, optional): Output directory. Default: `{protein_name}_results/`
- `pdb_id` (str, optional): PDB ID for structure mapping
- `max_sequences` (int): Maximum sequences to retrieve. Default: 50
- `min_conservation` (float): Conservation threshold for motifs. Range: [0, 1]. Default: 0.70
- `window_sizes` (list of int): Window sizes for scanning. Default: [7, 9, 11, 13, 15]
- `threads` (int): CPU threads for MAFFT. Default: 4
- `verbose` (bool): Print progress messages. Default: False

**Returns**:
- `AnalysisResults` object with attributes:
  - `.motifs`: List of dicts with motif information
  - `.conservation_scores`: NumPy array of per-position scores
  - `.alignment`: BioPython MultipleSeqAlignment
  - `.tree`: ETE3 Tree object
  - `.data`: Dict with summary statistics

**Methods on returned object**:
- `.summary()`: Print formatted summary
- `.export_json(path)`: Save results as JSON
- `.get_file(name)`: Get path to specific output file

**Example**:
```python
results = evomotif.analyze_protein("p53", "user@email.com")

# Access motifs
for motif in results.motifs:
    print(f"{motif['sequence']}: conservation={motif['conservation']:.3f}")

# Get files
print(f"Alignment: {results.get_file('alignment')}")
print(f"Tree: {results.get_file('tree')}")

# Export
results.export_json("p53_results.json")
```

### Module-Level Classes

#### `SequenceRetriever`

Retrieve sequences from NCBI.

```python
from evomotif.retrieval import SequenceRetriever

retriever = SequenceRetriever(email="your@email.com")
sequences = retriever.fetch_sequences(
    query="hemoglobin",
    max_results=100,
    database="protein"
)
```

**Methods**:
- `fetch_sequences(query, max_results=50, database="protein")`: Returns list of SeqRecord
- `filter_quality(sequences)`: Remove low-quality sequences
- `remove_redundancy(sequences, threshold=0.95)`: Cluster and select representatives

#### `SequenceAligner`

Perform multiple sequence alignment.

```python
from evomotif.alignment import SequenceAligner

aligner = SequenceAligner(threads=8)
alignment = aligner.align(
    sequences,
    output="aligned.fasta",
    method="mafft"
)
```

**Methods**:
- `align(sequences, output=None)`: Run MAFFT, returns MultipleSeqAlignment
- `check_quality(alignment)`: Calculate quality metrics

#### `ConservationScorer`

Calculate conservation scores.

```python
from evomotif.conservation import ConservationScorer

scorer = ConservationScorer()
conservation = scorer.calculate_conservation_scores(
    alignment,
    method="combined",
    weights=(0.5, 0.5)
)
```

**Methods**:
- `calculate_conservation_scores(alignment, method='combined', weights=(0.5, 0.5))`: Returns NumPy array
- `shannon_entropy(alignment)`: Shannon entropy scores
- `blosum_score(alignment)`: BLOSUM62 scores

**Options for `method`**:
- `"shannon"`: Shannon entropy only
- `"blosum"`: BLOSUM62 only
- `"combined"`: Weighted combination (recommended)

#### `MotifFinder`

Discover conserved motifs.

```python
from evomotif.motif_discovery import MotifFinder

finder = MotifFinder()
motifs = finder.find_motifs(
    conservation_scores,
    threshold=0.70,
    window_sizes=[7, 9, 11],
    max_overlap=0.5
)
```

**Methods**:
- `find_motifs(conservation, threshold, window_sizes)`: Returns list of motif dicts
- `validate_motifs(motifs, conservation)`: Apply statistical tests

#### `PhylogenyBuilder`

Build phylogenetic tree.

```python
from evomotif.phylogeny import PhylogenyBuilder

builder = PhylogenyBuilder()
tree = builder.build_tree(
    alignment,
    method="fasttree",
    output="tree.nwk"
)
```

**Methods**:
- `build_tree(alignment, method='fasttree')`: Returns ETE3 Tree
- `plot_tree(tree, output='tree.png')`: Visualize tree

#### `StructureMapper`

Map conservation to 3D structure.

```python
from evomotif.structure import StructureMapper

mapper = StructureMapper()
mapper.map_conservation_to_structure(
    pdb_file="1TUP.pdb",
    conservation_scores=conservation,
    output="conserved.pdb"
)
```

**Methods**:
- `map_conservation_to_structure(pdb_file, conservation, output)`: Write PDB with conservation in B-factor
- `get_alphafold_confidence(pdb_file, chain='A')`: Extract pLDDT scores from AlphaFold structure

---

## Implementation Guide

### Basic Workflow

```python
import evomotif

# One-line analysis
results = evomotif.analyze_protein("myprotein", "email@example.com")

# View results
print(results.summary())

# Export
results.export_json("results.json")
```

### Custom Workflow

```python
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifFinder
from evomotif.stats import StatisticalValidator
from evomotif.phylogeny import PhylogenyBuilder

# Step 1: Retrieve
retriever = SequenceRetriever(email="email@example.com")
sequences = retriever.fetch_sequences("p53", max_results=100)
print(f"Retrieved {len(sequences)} sequences")

# Step 2: Align
aligner = SequenceAligner(threads=8)
alignment = aligner.align(sequences, output="p53_aligned.fasta")
print(f"Alignment: {alignment.get_alignment_length()} positions")

# Step 3: Conservation
scorer = ConservationScorer()
conservation = scorer.calculate_conservation_scores(
    alignment,
    method="combined",
    weights=(0.6, 0.4)  # Emphasize Shannon
)

# Step 4: Find motifs
finder = MotifFinder()
motifs = finder.find_motifs(
    conservation,
    threshold=0.75,
    window_sizes=[7, 9, 11]
)
print(f"Found {len(motifs)} candidate motifs")

# Step 5: Validate
validator = StatisticalValidator()
validated_motifs = validator.validate_motifs(
    motifs,
    conservation,
    method="permutation",
    n_permutations=10000,
    fdr_alpha=0.05
)
print(f"{len(validated_motifs)} motifs pass statistical tests")

# Step 6: Phylogeny
builder = PhylogenyBuilder()
tree = builder.build_tree(alignment, output="p53_tree.nwk")

# Step 7: Results
for motif in validated_motifs:
    print(f"\nMotif: {motif['sequence']}")
    print(f"  Position: {motif['start']}-{motif['end']}")
    print(f"  Conservation: {motif['conservation']:.3f}")
    print(f"  P-value: {motif['p_value']:.2e}")
    print(f"  Effect size: {motif['effect_size']:.2f}")
```

### Batch Analysis

```python
import evomotif
import json

proteins = ["p53", "BRCA1", "EGFR", "MYC", "KRAS"]
results_summary = {}

for protein in proteins:
    print(f"\nAnalyzing {protein}...")
    
    results = evomotif.analyze_protein(
        protein,
        "email@example.com",
        output_dir=f"./results/{protein}"
    )
    
    results_summary[protein] = {
        'n_sequences': len(results.alignment),
        'n_motifs': len(results.motifs),
        'mean_conservation': float(results.data['mean_conservation']),
        'max_conservation': float(results.data['max_conservation'])
    }

# Save summary
with open("batch_summary.json", "w") as f:
    json.dump(results_summary, f, indent=2)

print("\nBatch analysis complete:")
for protein, data in results_summary.items():
    print(f"{protein}: {data['n_motifs']} motifs, "
          f"mean conservation={data['mean_conservation']:.3f}")
```

### Integration with AlphaFold

```python
from evomotif.structure import StructureMapper
import evomotif

# Analyze protein
results = evomotif.analyze_protein("p53", "email@example.com")

# Get AlphaFold structure (download from AlphaFold DB)
alphafold_pdb = "AF-P04637-F1-model_v4.pdb"

# Map conservation
mapper = StructureMapper()
mapper.map_conservation_to_structure(
    pdb_file=alphafold_pdb,
    conservation_scores=results.conservation_scores,
    output="p53_conserved.pdb"
)

# Extract pLDDT confidence
confidence = mapper.get_alphafold_confidence(alphafold_pdb, chain='A')

# Correlate conservation with confidence
import numpy as np
correlation = np.corrcoef(results.conservation_scores, confidence)[0, 1]
print(f"Conservation-confidence correlation: {correlation:.3f}")

# High correlation suggests conserved regions are structurally confident
```

---

## Use Cases

### 1. Mutagenesis Planning

**Scenario**: Design mutations to probe protein function

**Approach**:
- Identify highly conserved residues (conservation > 0.85)
- These are likely functionally critical
- Mutations here will likely disrupt function
- Conversely, low conservation residues (< 0.4) are safe mutation targets

**Example**:
```python
results = evomotif.analyze_protein("kinase_of_interest", "email@example.com")

# Find critical residues
critical_positions = []
for i, score in enumerate(results.conservation_scores, start=1):
    if score > 0.85:
        residue = results.alignment[0][i-1]  # Reference sequence
        critical_positions.append((i, residue, score))

print("Critical residues (avoid mutating):")
for pos, res, score in critical_positions:
    print(f"  Position {pos}: {res} (conservation={score:.3f})")

# Find variable residues (safe to mutate)
variable_positions = []
for i, score in enumerate(results.conservation_scores, start=1):
    if score < 0.4:
        residue = results.alignment[0][i-1]
        variable_positions.append((i, residue, score))

print("\nVariable residues (safe mutation targets):")
for pos, res, score in variable_positions:
    print(f"  Position {pos}: {res} (conservation={score:.3f})")
```

### 2. Disease Variant Interpretation

**Scenario**: Assess pathogenicity of missense mutations

**Approach**:
- Mutations in highly conserved positions more likely pathogenic
- Check if mutation alters physicochemical properties
- Compare to BLOSUM62 penalty

**Example**:
```python
results = evomotif.analyze_protein("BRCA1", "email@example.com")

# Analyze specific mutation: R175H (Arg -> His at position 175)
position = 175
conservation = results.conservation_scores[position - 1]
alignment_col = results.alignment[:, position - 1]

print(f"Position {position}:")
print(f"  Conservation: {conservation:.3f}")
print(f"  Observed residues: {set(alignment_col)}")

# Check BLOSUM62 penalty for R -> H
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")
penalty = blosum62[('R', 'H')]

print(f"  BLOSUM62(R->H): {penalty}")

# Interpretation
if conservation > 0.8:
    print("  WARNING: Highly conserved position, mutation likely pathogenic")
elif conservation > 0.6:
    print("  CAUTION: Moderately conserved, mutation may affect function")
else:
    print("  INFO: Variable position, mutation may be tolerated")
```

### 3. Functional Domain Annotation

**Scenario**: Identify functional domains in unannotated protein

**Approach**:
- Conserved motifs often correspond to functional domains
- Map motifs to known domain databases (Pfam, SMART)
- Use motif patterns to predict function

**Example**:
```python
results = evomotif.analyze_protein("unknown_protein", "email@example.com")

print(f"Discovered {len(results.motifs)} conserved motifs:")
for i, motif in enumerate(results.motifs, 1):
    print(f"\nMotif {i}:")
    print(f"  Sequence: {motif['sequence']}")
    print(f"  Position: {motif['start']}-{motif['end']}")
    print(f"  Conservation: {motif['conservation']:.3f}")
    print(f"  Length: {len(motif['sequence'])} residues")
    
    # Check for known motifs (manual or via pattern matching)
    sequence = motif['sequence']
    if 'DFG' in sequence:
        print("  -> Possible kinase catalytic motif")
    elif 'CXXC' in sequence:
        print("  -> Possible zinc-binding motif")
    elif sequence.count('C') >= 4:
        print("  -> Possible cysteine-rich domain")
```

### 4. Protein Engineering

**Scenario**: Design minimal functional construct

**Approach**:
- Identify conserved core regions
- Remove variable flanking regions
- Create truncated but functional variant

**Example**:
```python
results = evomotif.analyze_protein("large_protein", "email@example.com")

# Find conserved core (longest consecutive high-conservation region)
core_start, core_end = None, None
current_start = None
max_length = 0

for i, score in enumerate(results.conservation_scores):
    if score > 0.65:  # Threshold for "core"
        if current_start is None:
            current_start = i
    else:
        if current_start is not None:
            length = i - current_start
            if length > max_length:
                max_length = length
                core_start = current_start
                core_end = i
            current_start = None

print(f"Conserved core: positions {core_start}-{core_end} ({max_length} residues)")
print(f"Original protein: {len(results.conservation_scores)} residues")
print(f"Minimal construct: {core_start}-{core_end}")
print(f"Reduction: {100 * (1 - max_length/len(results.conservation_scores)):.1f}%")

# Extract core sequence
core_sequence = str(results.alignment[0][core_start:core_end].seq)
print(f"\nCore sequence:\n{core_sequence}")
```

### 5. Phylogenetic Analysis

**Scenario**: Understand evolutionary relationships

**Approach**:
- Use phylogenetic tree from EvoMotif
- Map conservation to tree branches
- Identify clade-specific conservation patterns

**Example**:
```python
results = evomotif.analyze_protein("cytochrome_c", "email@example.com")

# Analyze tree
tree = results.tree
print(f"Tree has {len(tree)} leaves")

# Find most distant sequences
leaf_names = [leaf.name for leaf in tree.get_leaves()]
distances = []
for i, leaf1 in enumerate(tree.get_leaves()):
    for leaf2 in tree.get_leaves()[i+1:]:
        dist = tree.get_distance(leaf1, leaf2)
        distances.append((leaf1.name, leaf2.name, dist))

distances.sort(key=lambda x: x[2], reverse=True)
print("\nMost divergent sequences:")
for name1, name2, dist in distances[:5]:
    print(f"  {name1} <-> {name2}: distance={dist:.3f}")
```

---

## Limitations

### 1. Alignment Quality Dependency

**Issue**: Conservation scores are only as good as the alignment

**Impact**:
- Misaligned regions produce spurious conservation patterns
- Gaps can artificially inflate/deflate conservation

**Mitigation**:
- Always inspect alignment quality metrics
- Consider manual curation for critical applications
- Use high-quality sequences (avoid partial or fragmented sequences)

**Check alignment quality**:
```python
from evomotif.alignment import AlignmentQuality

quality = AlignmentQuality(alignment)
print(f"Gap %: {quality.gap_percentage():.1f}")  # Should be < 20%
print(f"Mean identity: {quality.mean_pairwise_identity():.1f}%")  # 30-70% typical

# Flag problematic regions
problematic = quality.identify_problematic_regions(gap_threshold=0.5)
if problematic:
    print(f"Warning: {len(problematic)} regions with >50% gaps")
```

### 2. Taxonomic Bias

**Issue**: NCBI database overrepresents certain clades (vertebrates, model organisms)

**Impact**:
- May miss lineage-specific conservation patterns
- Bias toward well-studied proteins

**Mitigation**:
- Manually curate sequence set if studying specific lineages
- Use taxonomic filters in NCBI query
- Consider UniProt as alternative database

### 3. Short Proteins

**Issue**: Proteins < 50 residues have insufficient signal

**Impact**:
- High false positive rate
- Statistical tests underpowered

**Mitigation**:
- Use stricter thresholds (min_conservation > 0.80)
- Require larger sample sizes (max_sequences > 100)
- Validate with structural data

### 4. Highly Divergent Families

**Issue**: Distant homologs (< 20% identity) difficult to align

**Impact**:
- Poor alignment quality
- Unreliable conservation scores

**Mitigation**:
- Consider structure-based alignment (PROMALS3D, MATT)
- Focus on well-conserved domains
- Use profile-based methods (PSI-BLAST) for retrieval

### 5. Computational Cost

**Issue**: Large proteins with many sequences are slow

**Performance**:
- Alignment: O(N² × L) where N = sequences, L = length
- Tree building: O(N² × L)
- Memory: O(N × L)

**Example**: BRCA1 (1863 residues) with 500 sequences:
- Alignment: ~45 minutes
- Memory: ~10 GB

**Mitigation**:
- Reduce sequences (50-100 usually sufficient)
- Use fewer threads if memory-limited
- Run on HPC cluster for large projects

---

## Troubleshooting

### NCBI Errors

**Error**: `HTTP Error 429: Too Many Requests`

**Cause**: NCBI rate limiting

**Solution**:
```python
from evomotif.retrieval import SequenceRetriever

# Add delay between requests
retriever = SequenceRetriever(email="email@example.com", delay=1.0)
sequences = retriever.fetch_sequences("protein", max_results=100)
```

### MAFFT Not Found

**Error**: `Command 'mafft' not found`

**Solution**:
```bash
# Verify installation
which mafft

# If not found, install:
conda install -c bioconda mafft

# Or specify path explicitly
export PATH="/path/to/mafft/bin:$PATH"
```

### Low Conservation Scores

**Issue**: All conservation scores are low (< 0.5)

**Possible causes**:
1. Protein family is highly divergent
2. Retrieved sequences are not true homologs
3. Alignment quality is poor

**Diagnosis**:
```python
# Check alignment identity
from evomotif.alignment import AlignmentQuality
quality = AlignmentQuality(alignment)
print(f"Mean pairwise identity: {quality.mean_pairwise_identity():.1f}%")

# If < 20%, sequences may be too divergent
```

**Solution**:
- Use more stringent retrieval (reduce max_sequences)
- Filter by E-value or identity threshold
- Consider profile-based search

### No Motifs Found

**Issue**: No motifs pass statistical tests

**Possible causes**:
1. Threshold too strict
2. Protein has no conserved regions
3. Insufficient sequences

**Solution**:
```python
# Try lower threshold
results = evomotif.analyze_protein(
    "protein",
    "email@example.com",
    min_conservation=0.60  # Lower from default 0.70
)

# Or increase sequences
results = evomotif.analyze_protein(
    "protein",
    "email@example.com",
    max_sequences=150  # Increase from default 50
)
```

### Memory Errors

**Error**: `MemoryError` or system slowdown

**Cause**: Large alignment (many sequences × long protein)

**Solution**:
```python
# Reduce sequences
results = evomotif.analyze_protein(
    "large_protein",
    "email@example.com",
    max_sequences=50  # Reduce from 100+
)

# Or process in chunks (for batch analysis)
proteins = ["p1", "p2", "p3", ...]
for protein in proteins:
    results = evomotif.analyze_protein(protein, email, max_sequences=30)
    # Process results immediately, then release memory
    del results
```

---

## Performance Optimization

### Parallelization

```python
# Use more threads for alignment
results = evomotif.analyze_protein(
    "protein",
    "email@example.com",
    threads=16  # Scales well up to ~8-16 cores
)
```

### Caching

```python
# Save intermediate results for re-analysis with different parameters
from evomotif.pipeline import EvoMotifPipeline

pipeline = EvoMotifPipeline(cache_dir="./cache")

# First run: retrieval and alignment cached
results1 = pipeline.analyze("protein", "email", min_conservation=0.70)

# Second run: reuses cached alignment, only re-runs motif finding
results2 = pipeline.analyze("protein", "email", min_conservation=0.75)
```

### Batch Processing

```python
# Parallel batch processing
from multiprocessing import Pool
import evomotif

def analyze_protein_wrapper(args):
    protein, email = args
    return evomotif.analyze_protein(protein, email)

proteins = ["p53", "BRCA1", "EGFR", "MYC"]
args = [(p, "email@example.com") for p in proteins]

# Run in parallel (careful with memory)
with Pool(4) as pool:
    results_list = pool.map(analyze_protein_wrapper, args)
```

---

## Contact & Support

- **Issues**: https://github.com/yourusername/evomotif/issues
- **Email**: taha@example.com
- **Documentation**: https://github.com/yourusername/evomotif/docs

For bug reports, please include:
- EvoMotif version
- Python version
- Operating system
- Full error traceback
- Minimal reproducible example

---

## References

1. Shannon, C.E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27(3), 379-423.

2. Henikoff, S. & Henikoff, J.G. (1992). Amino acid substitution matrices from protein blocks. *PNAS*, 89(22), 10915-10919.

3. Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289-300.

4. Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772-780.

5. Price, M.N., Dehal, P.S., & Arkin, A.P. (2010). FastTree 2 – approximately maximum-likelihood trees for large alignments. *PLoS ONE*, 5(3), e9490.

6. Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*. 2nd ed. Hillsdale, NJ: Lawrence Erlbaum.

7. Landau, M. et al. (2005). ConSurf 2005: the projection of evolutionary conservation scores of residues on protein structures. *Nucleic Acids Research*, 33, W299-W302.
