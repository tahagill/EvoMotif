# EvoMotif: Complete Technical Documentation

**Version:** 1.0.0  
**Last Updated:** December 23, 2025

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation & Setup](#installation--setup)
3. [Quick Start](#quick-start)
4. [Core Algorithms & Logic](#core-algorithms--logic)
5. [Statistical Methods](#statistical-methods)
6. [Module Reference](#module-reference)
7. [Output Formats](#output-formats)
8. [Advanced Usage](#advanced-usage)
9. [Best Practices](#best-practices)
10. [Troubleshooting](#troubleshooting)

---

## 1. Introduction

### What is EvoMotif?

EvoMotif is a comprehensive Python library for discovering evolutionarily conserved protein motifs through multi-species sequence analysis. It automates the entire workflow from sequence retrieval to statistical validation and 3D structural mapping.

### Key Features

- **One-Line Analysis**: Complete protein analysis with a single function call
- **Scientifically Rigorous**: Statistical validation with multiple testing correction
- **Publication Ready**: Generates figures, tables, and supplementary data
- **AlphaFold Integration**: Correlate conservation with structural confidence
- **Flexible**: Use high-level API or individual modules for custom workflows

### Scientific Rationale

**Why study conserved motifs?**

Evolutionarily conserved regions in proteins indicate:
- **Functional importance**: Essential for protein activity
- **Structural constraints**: Required for proper folding
- **Interaction sites**: Binding domains for other molecules
- **Disease relevance**: Mutations in conserved regions often pathogenic

**EvoMotif's Approach:**

1. **Multi-species comparison** reduces noise and increases signal
2. **Combined metrics** (Shannon entropy + BLOSUM62) capture both variability and functional constraints
3. **Statistical validation** ensures discoveries are not random patterns
4. **3D mapping** connects sequence conservation to structural context

---

## 2. Installation & Setup

### System Requirements

- **OS**: Linux, macOS, or Windows (WSL recommended)
- **Python**: 3.8 or higher
- **RAM**: 8GB minimum, 16GB recommended
- **Storage**: 2GB for software + space for results

### External Dependencies

EvoMotif requires these bioinformatics tools:

1. **MAFFT** (Multiple sequence alignment)
   ```bash
   # Ubuntu/Debian
   sudo apt-get install mafft
   
   # macOS
   brew install mafft
   
   # Conda
   conda install -c bioconda mafft
   ```

2. **FastTree** (Phylogenetic tree inference)
   ```bash
   # Ubuntu/Debian
   sudo apt-get install fasttree
   
   # macOS
   brew install fasttree
   
   # Conda
   conda install -c bioconda fasttree
   ```

3. **DSSP** (Optional - for secondary structure)
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
git clone https://github.com/yourusername/EvoMotif.git
cd EvoMotif

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install EvoMotif
pip install -e .

# Verify installation
python -c "import evomotif; print('EvoMotif installed successfully!')"
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

---

## 3. Quick Start

### Basic Analysis (2 Lines)

```python
import evomotif

# Analyze any protein by name
results = evomotif.analyze_protein("hemoglobin", "your@email.com")

# View summary
print(results.summary())
```

**Output:**
```
============================================================
EvoMotif Analysis Results: hemoglobin
============================================================
Sequences analyzed: 24
Consecutive motifs found: 9
Conserved positions: 73
Mean conservation: 0.642
Max conservation: 1.000

Top 5 motifs:
  1. GAEAL (pos 26-30, cons=0.735)
  2. HGKKV (pos 43-47, cons=0.843)
  3. SDLHA (pos 51-55, cons=0.688)
  ...

📁 Results saved to: hemoglobin_results/
============================================================
```

### With 3D Structure

```python
# Map conservation to PDB structure
results = evomotif.analyze_protein(
    protein_name="p53",
    email="your@email.com",
    pdb_id="1TUP",  # DNA-binding domain
    min_conservation=0.70
)

# Access structure file
structure_file = results.get_file('structure')
```

### Custom Parameters

```python
results = evomotif.analyze_protein(
    protein_name="BRCA1",
    email="your@email.com",
    output_dir="./my_analysis",
    max_sequences=100,        # More sequences = better statistics
    min_conservation=0.65,    # Lower threshold = more motifs
    threads=8,                # Parallel processing
    verbose=True             # Show detailed progress
)
```

### Working with Results

```python
# Access data programmatically
motifs = results.motifs
conserved_pos = results.conserved_positions
conservation = results.conservation_scores

# Export to JSON
results.export_json("my_results.json")

# Get specific files
fasta_file = results.get_file('sequences')
alignment_file = results.get_file('alignment')
tree_file = results.get_file('tree')

# Iterate over motifs
for motif in results.motifs:
    print(f"Motif at {motif['start']}-{motif['end']}: {motif['sequence']}")
    print(f"  Conservation: {motif['conservation']:.3f}")
    print(f"  Sequences: {len(motif['sequences'])}")
```

---

## 4. Core Algorithms & Logic

### 4.1 Sequence Retrieval Strategy

**Algorithm:** Intelligent NCBI querying with redundancy filtering

**Logic:**
1. Query NCBI Protein database with protein name
2. Filter by taxonomy (prioritize model organisms)
3. Remove redundant sequences (>95% identity)
4. Ensure diverse taxonomic sampling

**Why This Approach?**
- **Diversity**: Captures evolutionary breadth
- **Quality**: Model organisms have well-annotated sequences
- **Statistical power**: Multiple sequences improve conservation detection
- **Efficiency**: Redundancy removal reduces computational cost

**Implementation:**
```python
# Retrieval module logic
def fetch_sequences(protein_name, email, max_sequences=50):
    """
    Smart retrieval with quality filters.
    
    Steps:
    1. Search NCBI with protein name
    2. Fetch sequences in batches
    3. Remove sequences with ambiguous residues (X, B, Z)
    4. Cluster by similarity (CD-HIT-like approach)
    5. Select representative sequences
    6. Return diverse, high-quality set
    """
```

**Parameters Explained:**
- `max_sequences`: More sequences improve statistics but increase runtime
  - Minimum: 10 (basic analysis)
  - Recommended: 50-100 (robust analysis)
  - Maximum: 500 (comprehensive, slow)

### 4.2 Conservation Scoring Algorithm

**Combined Metric Approach**

EvoMotif uses **two complementary metrics** to measure conservation:

#### A. Shannon Entropy (Information Theory)

**What It Measures:** Variability/uncertainty at each position

**Mathematical Formula:**
$$H(i) = -\sum_{a \in \text{amino acids}} p_a(i) \log_2 p_a(i)$$

**Normalized Conservation Score:**
$$C_{\text{shannon}}(i) = 1 - \frac{H(i)}{\log_2(20)}$$

**Interpretation:**
- **Score = 1.0**: All sequences have identical amino acid (perfect conservation)
- **Score = 0.5**: Moderate variability
- **Score = 0.0**: Maximum variability (all 20 amino acids equally likely)

**Why Shannon Entropy?**
- Quantifies information content
- Sensitive to rare variants
- Mathematical foundation in information theory
- No biological assumptions

**Example:**
```
Position 45: [H, H, H, H, H, H, H, H, H, R]
- H appears 90%, R appears 10%
- H(45) = -(0.9×log₂(0.9) + 0.1×log₂(0.1)) = 0.469 bits
- C_shannon(45) = 1 - 0.469/4.322 = 0.891 (highly conserved)
```

#### B. BLOSUM62 Score (Evolutionary Perspective)

**What It Measures:** Evolutionary substitution patterns

**Mathematical Formula:**
$$B(i) = \frac{1}{n(n-1)} \sum_{j<k} \text{BLOSUM62}(a_j, a_k)$$

Where:
- $a_j, a_k$ = amino acids at position $i$ in sequences $j, k$
- BLOSUM62 = empirical substitution matrix from aligned protein blocks

**Normalized Score:**
$$B_{\text{norm}}(i) = \frac{B(i) - B_{\text{min}}}{B_{\text{max}} - B_{\text{min}}}$$

**BLOSUM62 Matrix Philosophy:**
- Positive scores: Functionally similar (conservative substitutions)
- Negative scores: Functionally dissimilar (radical substitutions)
- Example: H ↔ R = -1 (similar: both basic), H ↔ D = -3 (different: basic vs acidic)

**Why BLOSUM62?**
- Captures functional constraints
- Based on real evolutionary data
- Distinguishes conservative vs. radical changes
- Widely validated in bioinformatics

**Example:**
```
Position 45: [H, H, H, H, H, H, H, H, H, R]
Pairwise scores:
- BLOSUM62(H, H) = 8 (identical)
- BLOSUM62(H, R) = 0 (neutral substitution)
Average score ≈ 7.8 (high, indicates conservative evolution)
B_norm(45) = 0.78 (functionally conserved)
```

#### C. Combined Conservation Score

**Formula:**
$$C_{\text{final}}(i) = 0.5 \times C_{\text{shannon}}(i) + 0.5 \times B_{\text{norm}}(i)$$

**Why Combine Both?**
1. **Shannon**: Detects identical residues (structural/catalytic importance)
2. **BLOSUM**: Detects functional equivalence (e.g., Leu ↔ Ile, both hydrophobic)
3. **Together**: Captures both strict and functional conservation

**Real Example (from hemoglobin analysis):**
```
Position 59 (Histidine - heme binding):
- Shannon score: 0.93 (mostly His, few Tyr)
- BLOSUM score: 0.85 (His/Tyr both aromatic, functionally similar)
- Combined: 0.89 (highly conserved + functionally important)
```

**Implementation:**
```python
def calculate_conservation(alignment):
    """
    Calculate combined conservation scores.
    
    Algorithm:
    1. For each alignment column:
       a. Calculate Shannon entropy
       b. Calculate average BLOSUM62 score
       c. Normalize both to [0, 1]
       d. Combine with equal weights
    2. Return final scores (0 = variable, 1 = conserved)
    """
    shannon_scores = calculate_shannon_entropy(alignment)
    blosum_scores = calculate_blosum_score(alignment)
    
    # Equal weighting (can be tuned for specific applications)
    combined = 0.5 * shannon_scores + 0.5 * blosum_scores
    
    return combined
```

### 4.3 Motif Discovery Algorithm

**Sliding Window Approach with Adaptive Thresholds**

**Algorithm Steps:**

1. **Sliding Window Scanning**
   ```
   Alignment: M-S-L-S-D-K-G-K-A-T-V-R-A-I-W-G-K-I-G
   Conservation: 0.8 0.7 0.8 0.7 0.8 0.7 0.9 0.7 0.8 0.9...
   
   Window size 5:
   [M-S-L-S-D] → avg = 0.76 → candidate
   [S-L-S-D-K] → avg = 0.74 → candidate
   [L-S-D-K-G] → avg = 0.78 → candidate
   ...
   ```

2. **Threshold Filtering**
   - Default threshold: 0.70 (tunable)
   - Only windows with avg conservation ≥ threshold pass

3. **Overlap Resolution**
   ```
   Overlapping candidates:
   [M-S-L-S-D] (0.76)
   [S-L-S-D-K] (0.74)
   [L-S-D-K-G] (0.78) ← Keep highest
   ```

4. **Consecutive Motif Merging**
   ```
   Adjacent high-scoring windows merge:
   [A-I-W-G-K] (0.82) + [W-G-K-I-G] (0.79) 
   → [A-I-W-G-K-I-G] (0.805)
   ```

5. **Statistical Validation**
   - Permutation test: Shuffle conservation scores
   - Calculate p-value: P(random ≥ observed)
   - FDR correction: Control for multiple testing

**Why This Algorithm?**

- **Sliding windows**: Detect local conservation patterns
- **Multiple sizes**: Capture motifs of varying length (7-15 residues typical)
- **Overlap resolution**: Avoids redundant motifs
- **Merging**: Finds extended conserved regions
- **Validation**: Ensures statistical significance

**Parameters:**

```python
# Motif discovery parameters
window_sizes = [7, 9, 11, 13, 15]  # Try multiple lengths
min_conservation = 0.70            # Threshold for candidates
max_std = 0.15                     # Maximum variability within motif
min_gap_free = 0.70                # Require ≥70% sequences without gaps
```

**Parameter Tuning Guide:**
- **Higher `min_conservation`** (0.75-0.80): Stricter, fewer motifs, higher confidence
- **Lower `min_conservation`** (0.60-0.65): Permissive, more motifs, lower confidence
- **Window sizes**: Smaller (5-7) for short motifs, larger (13-15) for domains

### 4.4 Statistical Validation Logic

**Multi-Level Validation Framework**

#### Level 1: Permutation Testing

**Null Hypothesis:** Observed conservation pattern is random

**Test Procedure:**
```python
def permutation_test(observed_conservation, n_permutations=10000):
    """
    Test if motif conservation is statistically significant.
    
    Steps:
    1. Calculate observed statistic (mean conservation of motif)
    2. Generate null distribution:
       - Randomly shuffle conservation scores
       - Extract same-length window
       - Calculate mean conservation
       - Repeat 10,000 times
    3. Calculate p-value: fraction of null ≥ observed
    4. Reject null if p < 0.05
    """
    observed_mean = np.mean(observed_conservation)
    
    null_distribution = []
    for _ in range(n_permutations):
        shuffled = np.random.permutation(all_conservation_scores)
        random_window = shuffled[:len(observed_conservation)]
        null_distribution.append(np.mean(random_window))
    
    p_value = np.mean(null_distribution >= observed_mean)
    return p_value
```

**Why Permutation Tests?**
- **Non-parametric**: No assumptions about data distribution
- **Exact**: p-values are exact given the data
- **Powerful**: Sensitive to patterns in real data

#### Level 2: Multiple Testing Correction

**Problem:** Testing many motifs increases false positives

**Solution:** False Discovery Rate (FDR) control via Benjamini-Hochberg

**Algorithm:**
```python
def fdr_correction(p_values, alpha=0.05):
    """
    Control false discovery rate across multiple tests.
    
    Benjamini-Hochberg procedure:
    1. Sort p-values: p₁ ≤ p₂ ≤ ... ≤ pₘ
    2. Find largest i where: pᵢ ≤ (i/m) × α
    3. Reject H₀ for all tests 1 to i
    
    Result: Expected FDR ≤ α (e.g., ≤5% false discoveries)
    """
```

**Why FDR vs. Bonferroni?**
- **Bonferroni**: Too conservative, misses true motifs
- **FDR**: Balances discovery and false positives
- **Appropriate**: When finding multiple true signals is expected

#### Level 3: Effect Size Calculation

**Cohen's d for Conservation Differences:**

$$d = \frac{\mu_{\text{motif}} - \mu_{\text{background}}}{\sigma_{\text{pooled}}}$$

**Interpretation:**
- |d| < 0.2: Small effect (weak conservation)
- |d| = 0.5: Medium effect (moderate conservation)
- |d| > 0.8: Large effect (strong conservation)

**Why Effect Sizes?**
- **Practical significance**: p-value alone doesn't indicate importance
- **Standardized**: Comparable across studies
- **Publication quality**: Recommended by statistical guidelines

---

## 5. Statistical Methods

### 5.1 Conservation Metrics Summary

| Metric | Formula | Range | Interpretation |
|--------|---------|-------|----------------|
| Shannon Entropy | $H = -\sum p \log_2 p$ | 0 (conserved) to 4.32 (variable) | Information content |
| Shannon Conservation | $C = 1 - H/H_{\max}$ | 0 (variable) to 1 (conserved) | Inverse entropy |
| BLOSUM62 Score | $\overline{\text{BLOSUM}(i,j)}$ | -4 (dissimilar) to 11 (identical) | Evolutionary constraint |
| Combined Score | $0.5C_s + 0.5B_n$ | 0 to 1 | Comprehensive conservation |
| Gap Frequency | $f_{\text{gap}} = n_{\text{gaps}}/n_{\text{seqs}}$ | 0 (no gaps) to 1 (all gaps) | Structural flexibility |

### 5.2 Hypothesis Tests Used

| Test | Purpose | When Used | Why Chosen |
|------|---------|-----------|------------|
| **Permutation Test** | Motif significance | Every discovered motif | Non-parametric, exact |
| **Mann-Whitney U** | Compare two groups | Variant enrichment | Robust to non-normality |
| **Fisher's Exact** | Categorical associations | Variant distribution | Exact for small samples |
| **Pearson Correlation** | Linear relationships | Conservation vs. structure | Standard for continuous data |
| **Kolmogorov-Smirnov** | Distribution comparison | QC checks | Sensitive to differences |

### 5.3 Multiple Testing Corrections

**Available Methods:**

1. **Bonferroni**: $\alpha_{\text{adj}} = \alpha / m$
   - Use: Small number of tests (<10)
   - Conservative: Controls family-wise error rate

2. **FDR (Benjamini-Hochberg)**: Adaptive threshold
   - Use: Many tests (10-1000s)
   - Recommended: Balances power and false positives

3. **FDR (Benjamini-Yekutieli)**: For dependent tests
   - Use: Tests are not independent
   - Conservative than B-H but valid for dependencies

**Implementation:**
```python
from statsmodels.stats.multitest import multipletests

reject, adjusted_p, _, _ = multipletests(
    p_values,
    alpha=0.05,
    method='fdr_bh'  # or 'bonferroni', 'fdr_by'
)
```

### 5.4 Quality Control Metrics

**Alignment Quality:**
```python
# Gap percentage (should be <20%)
gap_pct = (n_gaps / (n_seqs * length)) * 100

# Sequence identity (should be 30-90%)
avg_identity = mean_pairwise_identity(alignment)

# Coverage (should be >80%)
coverage = n_aligned_positions / query_length
```

**Conservation Quality:**
```python
# Dynamic range (should span 0.3-1.0)
conservation_range = max(scores) - min(scores)

# Distribution shape (should be bimodal)
hist, bins = np.histogram(scores, bins=20)

# Correlation with structure (if available)
r, p = pearsonr(conservation, bfactor)
```

---

## 6. Module Reference

### 6.1 Pipeline Module (`evomotif.pipeline`)

**Main Entry Point:**

```python
from evomotif import analyze_protein

results = analyze_protein(
    protein_name: str,              # Protein name or identifier
    email: str,                     # NCBI requires email
    output_dir: str = None,         # Output directory (default: protein_name/)
    pdb_id: str = None,             # PDB ID for structure mapping
    max_sequences: int = 50,        # Number of sequences to retrieve
    min_conservation: float = 0.70, # Conservation threshold
    threads: int = 4,               # CPU cores to use
    verbose: bool = False,          # Print detailed progress
    check_deps: bool = True         # Check dependencies first
) -> AnalysisResults
```

**AnalysisResults Object:**

```python
class AnalysisResults:
    # Attributes
    .protein: str                    # Protein name
    .n_sequences: int                # Number of sequences analyzed
    .motifs: List[Dict]              # Discovered motifs
    .conserved_positions: List[Dict] # Highly conserved positions
    .conservation_scores: ndarray    # Conservation per position
    .output_dir: Path                # Results directory
    .files: Dict[str, Path]          # Output file paths
    
    # Methods
    .summary() -> str                # Human-readable summary
    .get_file(type) -> Path          # Get specific file path
    .export_json(path)               # Export to JSON
```

### 6.2 Retrieval Module (`evomotif.retrieval`)

**SequenceRetriever Class:**

```python
from evomotif.retrieval import SequenceRetriever

retriever = SequenceRetriever(
    email="your@email.com",          # Required by NCBI
    api_key=None,                    # Optional: faster queries
    max_retries=3,                   # Retry failed requests
    delay=0.5                        # Delay between requests
)

# Fetch sequences
sequences = retriever.fetch_sequences(
    query="hemoglobin alpha",
    max_results=50,
    database="protein"
)

# Save to FASTA
retriever.save_fasta(sequences, "output.fasta")
```

### 6.3 Alignment Module (`evomotif.alignment`)

**SequenceAligner Class:**

```python
from evomotif.alignment import SequenceAligner

aligner = SequenceAligner(
    tool="mafft",                    # Alignment tool
    threads=4,                       # CPU cores
    algorithm="auto"                 # MAFFT algorithm
)

# Align sequences
alignment = aligner.align(
    sequences="input.fasta",
    output="aligned.fasta",
    remove_gaps=False
)

# Alignment statistics
stats = aligner.get_alignment_stats(alignment)
# Returns: n_sequences, length, gap_percentage, identity
```

### 6.4 Conservation Module (`evomotif.conservation`)

**ConservationScorer Class:**

```python
from evomotif.conservation import ConservationScorer

scorer = ConservationScorer()

# Calculate conservation
conservation = scorer.calculate_conservation_scores(
    alignment,
    method="combined",               # "shannon", "blosum62", or "combined"
    weights=(0.5, 0.5),             # (shannon_weight, blosum_weight)
    ignore_gaps=True
)

# Identify conserved positions
conserved = scorer.find_conserved_positions(
    alignment,
    conservation,
    threshold=0.70,
    min_gap_free=0.70
)
```

### 6.5 Motif Discovery Module (`evomotif.motif_discovery`)

**MotifDiscoverer Class:**

```python
from evomotif.motif_discovery import MotifDiscoverer

discoverer = MotifDiscoverer(
    window_sizes=[7, 9, 11, 13, 15],
    min_conservation=0.70,
    max_std=0.15,
    min_gap_free=0.70
)

# Discover consecutive motifs
motifs = discoverer.find_consecutive_motifs(
    alignment,
    conservation_scores
)

# Find scattered conserved residues
scattered = discoverer.find_scattered_residues(
    alignment,
    conservation_scores,
    min_distance=3,
    max_distance=20
)
```

### 6.6 Statistical Analysis Module (`evomotif.stats`)

**StatisticalAnalyzer Class:**

```python
from evomotif.stats import StatisticalAnalyzer

analyzer = StatisticalAnalyzer(random_seed=42)

# Permutation test
p_value = analyzer.motif_permutation_test(
    motif_conservation,
    all_conservation_scores,
    n_permutations=10000
)

# Multiple testing correction
reject, adj_p = analyzer.multiple_testing_correction(
    p_values,
    method='fdr_bh',
    alpha=0.05
)

# Effect size
effect = analyzer.calculate_effect_size(
    group1_scores,
    group2_scores,
    method='cohens_d'
)
```

### 6.7 Phylogenetic Module (`evomotif.phylogeny`)

**PhylogeneticAnalyzer Class:**

```python
from evomotif.phylogeny import PhylogeneticAnalyzer

phylo = PhylogeneticAnalyzer(
    method="ml",                     # Maximum likelihood
    bootstrap=100                    # Bootstrap replicates
)

# Build tree
tree = phylo.build_tree(
    alignment_file="aligned.fasta",
    output_file="tree.nwk"
)

# Visualize tree
phylo.plot_tree(
    tree_file="tree.nwk",
    output_file="tree.png",
    show_labels=True
)
```

### 6.8 Structure Module (`evomotif.structure`)

**StructureMapper Class:**

```python
from evomotif.structure import StructureMapper

mapper = StructureMapper()

# Download PDB structure
structure = mapper.download_structure(
    pdb_id="1TUP",
    output_dir="structures/"
)

# Map conservation to structure
mapper.map_conservation_to_structure(
    structure=structure,
    conservation_scores=conservation,
    sequence_alignment=alignment,
    output_file="conserved_structure.pdb"
)

# Extract AlphaFold confidence
confidence = mapper.get_alphafold_confidence(
    structure=alphafold_structure,
    chain_id='A'
)
```

---

## 7. Output Formats

### 7.1 File Structure

```
protein_name_results/
├── protein_name_sequences.fasta        # Retrieved sequences
├── protein_name_aligned.fasta          # Multiple sequence alignment
├── protein_name_conservation.json      # Conservation scores
├── protein_name_tree.nwk              # Phylogenetic tree (Newick)
├── protein_name_tree.png              # Tree visualization
├── conserved_positions.json            # High-conservation positions
├── protein_name_summary.json          # Complete analysis summary
├── motifs/
│   ├── motif_1.fasta                  # Individual motif sequences
│   ├── motif_2.fasta
│   └── ...
└── structure/
    └── conserved_structure.pdb        # PDB with conservation in B-factor
```

### 7.2 JSON Formats

**Conservation JSON:**
```json
{
  "protein": "hemoglobin alpha",
  "alignment_length": 143,
  "conservation_scores": [0.799, 0.569, 0.766, ...],
  "mean_conservation": 0.642,
  "median_conservation": 0.637,
  "std_conservation": 0.139
}
```

**Conserved Positions JSON:**
```json
[
  {
    "position": 14,
    "conservation": 0.999,
    "gap_frequency": 0.0,
    "amino_acid": "W",
    "annotation": "Highly conserved tryptophan"
  },
  ...
]
```

**Summary JSON:**
```json
{
  "protein": "hemoglobin alpha",
  "n_sequences": 24,
  "alignment_stats": {
    "alignment_length": 143,
    "gap_percentage": 0.44,
    "mean_identity": 0.73
  },
  "motifs": [
    {
      "motif_id": 1,
      "start": 26,
      "end": 30,
      "sequence": "GAEAL",
      "conservation": 0.735,
      "p_value": 0.0001,
      "significant": true
    }
  ],
  "statistics": {
    "n_motifs": 9,
    "n_conserved_positions": 73,
    "mean_conservation": 0.642
  }
}
```

### 7.3 FASTA Formats

**Sequences FASTA:**
```
>NP_001299611.1 hemoglobin subunit alpha [Ictidomys tridecemlineatus]
MVLSPADKNNVKACWEKIGGHGAAYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVQGHG
KKVADALANAAAHVDDLPGALSTLSDLHAHKLRVDPVNFKLLSHCLLVTLAAHHPAEFTP
AVHASLDKFLASVSTVLTSKYC
```

**Aligned FASTA:**
```
>seq1
MVLSPADKNNVKACWEKIGGHGAAYGAEAL---ERMFLSFPTTKTYFPHFDLSHGSAQVQGHG
>seq2
MSLSDTDKAVVRALWGKISSRSDDIGAEAL---GRMLTVYPQTKTYFSHWADLSPGSAPVKKH
```

**Motif FASTA:**
```
>motif1_seq1 pos_26-30_cons_0.735
GAEAL
>motif1_seq2 pos_26-30_cons_0.735
GAEAL
>motif1_seq3 pos_26-30_cons_0.735
GNNAL
```

---

## 8. Advanced Usage

### 8.1 Custom Conservation Metrics

```python
from evomotif.conservation import ConservationScorer

scorer = ConservationScorer()

# Shannon entropy only
shannon = scorer.calculate_shannon_entropy(alignment)

# BLOSUM62 only
blosum = scorer.calculate_blosum_score(alignment)

# Custom weighting
custom = scorer.calculate_conservation_scores(
    alignment,
    method="combined",
    weights=(0.7, 0.3)  # 70% Shannon, 30% BLOSUM
)
```

### 8.2 Advanced Motif Discovery

```python
from evomotif.motif_discovery import MotifDiscoverer

# Strict criteria for high-confidence motifs
strict_discoverer = MotifDiscoverer(
    window_sizes=[9, 11],           # Focus on medium-length
    min_conservation=0.80,          # Higher threshold
    max_std=0.10,                   # Lower variability
    min_gap_free=0.90               # Fewer gaps allowed
)

# Permissive criteria for exploratory analysis
permissive_discoverer = MotifDiscoverer(
    window_sizes=[5, 7, 9, 11, 13, 15, 17, 19],
    min_conservation=0.60,
    max_std=0.20,
    min_gap_free=0.60
)
```

### 8.3 Custom Statistical Tests

```python
from evomotif.stats import StatisticalAnalyzer

analyzer = StatisticalAnalyzer()

# Custom permutation test
def my_statistic(data):
    return np.percentile(data, 90)  # 90th percentile

result = analyzer.permutation_test(
    observed_statistic=0.85,
    data=conservation_scores,
    statistic_func=my_statistic,
    n_permutations=10000,
    alternative='greater'
)
```

### 8.4 Batch Analysis

```python
proteins = ["p53", "BRCA1", "EGFR", "TP53", "MYC"]
all_results = {}

for protein in proteins:
    try:
        results = evomotif.analyze_protein(
            protein,
            "your@email.com",
            output_dir=f"batch_results/{protein}"
        )
        all_results[protein] = results
        print(f"✓ {protein}: {len(results.motifs)} motifs found")
    except Exception as e:
        print(f"✗ {protein}: {str(e)}")

# Compare results
for protein, results in all_results.items():
    print(f"{protein}: {results.data['mean_conservation']:.3f} mean conservation")
```

### 8.5 Integration with AlphaFold

```python
from evomotif.structure import StructureMapper

mapper = StructureMapper()

# Load AlphaFold structure
af_structure = mapper.download_structure(
    pdb_id="AF-P04637",  # AlphaFold ID for p53
    output_dir="alphafold_structures/"
)

# Extract pLDDT confidence scores
confidence = mapper.get_alphafold_confidence(af_structure, chain_id='A')

# Correlate with conservation
from scipy.stats import pearsonr

# Align positions
common_positions = set(confidence.keys()) & set(range(len(conservation)))
conf_values = [confidence[p] for p in sorted(common_positions)]
cons_values = [conservation[p] for p in sorted(common_positions)]

r, p_value = pearsonr(conf_values, cons_values)
print(f"Conservation vs. pLDDT: r={r:.3f}, p={p_value:.2e}")
```

---

## 9. Best Practices

### 9.1 Sequence Selection

**Recommendations:**
- **Minimum 10 sequences**: Basic statistics
- **Optimal 50-100 sequences**: Robust conservation detection
- **Maximum 500 sequences**: Comprehensive but slow

**Taxonomic Diversity:**
- Include sequences from multiple domains/kingdoms
- Avoid over-representation of closely related species
- Model organisms provide high-quality annotations

### 9.2 Parameter Tuning

**Conservation Threshold:**
```python
# Strict (high confidence, fewer motifs)
min_conservation = 0.80

# Standard (balanced)
min_conservation = 0.70

# Permissive (exploratory, more motifs)
min_conservation = 0.60
```

**Motif Length:**
```python
# Short functional sites (kinase sites, etc.)
window_sizes = [5, 7, 9]

# Medium motifs (binding sites)
window_sizes = [7, 9, 11, 13, 15]

# Long domains
window_sizes = [13, 15, 17, 19, 21]
```

### 9.3 Quality Control

**Check Alignment Quality:**
```python
# Gap percentage should be <20%
if stats['gap_percentage'] > 20:
    print("WARNING: High gap percentage, consider:")
    print("- Using fewer sequences")
    print("- Filtering for length")
    print("- Different alignment parameters")

# Sequence identity should be 30-90%
if stats['mean_identity'] < 30:
    print("WARNING: Low sequence identity, may not be homologous")
if stats['mean_identity'] > 90:
    print("WARNING: High sequence identity, insufficient diversity")
```

**Validate Conservation Scores:**
```python
# Check dynamic range
conservation_range = max(conservation) - min(conservation)
if conservation_range < 0.3:
    print("WARNING: Low conservation range, all positions similar")

# Check for bimodal distribution (expected)
import matplotlib.pyplot as plt
plt.hist(conservation, bins=20)
plt.xlabel('Conservation Score')
plt.ylabel('Frequency')
plt.title('Conservation Distribution')
plt.show()
```

### 9.4 Interpretation Guidelines

**Conservation Score Interpretation:**
- **0.90-1.00**: Catalytic residues, structural cores, critical interactions
- **0.75-0.90**: Functional sites, binding pockets, domain interfaces
- **0.60-0.75**: Moderately conserved, potential functional importance
- **0.40-0.60**: Variable regions, loops, linkers
- **0.00-0.40**: Highly variable, likely non-functional

**Statistical Significance:**
- **p < 0.001**: Very strong evidence
- **p < 0.01**: Strong evidence
- **p < 0.05**: Moderate evidence
- **p ≥ 0.05**: Insufficient evidence

**Effect Sizes:**
- **|d| < 0.2**: Small effect, limited practical importance
- **|d| = 0.5**: Medium effect, notable difference
- **|d| > 0.8**: Large effect, substantial conservation

---

## 10. Troubleshooting

### 10.1 Common Issues

**Issue: No sequences retrieved**
```
Solution:
1. Check protein name spelling
2. Try alternative names (gene symbol vs. protein name)
3. Verify NCBI email is valid
4. Check internet connection
5. Try with organism filter: "hemoglobin alpha AND human[organism]"
```

**Issue: Alignment fails**
```
Solution:
1. Verify MAFFT is installed: `mafft --version`
2. Check sequence file is not empty
3. Ensure sequences are protein (not DNA)
4. Try fewer sequences if memory error
5. Check file permissions
```

**Issue: No motifs found**
```
Solution:
1. Lower min_conservation threshold (try 0.60)
2. Use more sequences (target 50-100)
3. Check conservation scores distribution
4. Try different window sizes
5. Verify alignment quality (gap percentage)
```

**Issue: All p-values are 1.0**
```
Solution:
1. Increase n_permutations (try 50,000)
2. Check conservation scores are not uniform
3. Verify motif discovery found real patterns
4. Review statistical test implementation
```

### 10.2 Performance Optimization

**Speed Up Analysis:**
```python
# Use more threads
results = evomotif.analyze_protein(
    protein,
    email,
    threads=8  # Use all available cores
)

# Reduce sequences for quick tests
results = evomotif.analyze_protein(
    protein,
    email,
    max_sequences=20  # Faster but less robust
)

# Skip optional analyses
pipeline = evomotif.EvoMotifPipeline()
results = pipeline.analyze(
    protein,
    email,
    skip_tree=True,      # Skip phylogenetic tree
    skip_structure=True  # Skip structure mapping
)
```

**Memory Optimization:**
```python
# For large analyses (>500 sequences)
import gc

for chunk in sequence_chunks:
    results = analyze_chunk(chunk)
    process_results(results)
    gc.collect()  # Free memory
```

### 10.3 Error Messages

| Error | Meaning | Solution |
|-------|---------|----------|
| `ModuleNotFoundError: evomotif` | Not installed | `pip install -e .` |
| `FileNotFoundError: mafft` | MAFFT not in PATH | Install MAFFT |
| `HTTPError: 429` | Too many NCBI requests | Add delay, get API key |
| `ValueError: Empty alignment` | No sequences aligned | Check input sequences |
| `MemoryError` | Insufficient RAM | Reduce sequence count |

### 10.4 Getting Help

**Check Documentation:**
1. Read error message carefully
2. Search this guide (Ctrl+F)
3. Review example code in `examples/`
4. Check STATISTICS_METHODS.md for math details

**Debug Mode:**
```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Now run analysis - will print detailed logs
results = evomotif.analyze_protein(protein, email, verbose=True)
```

**Report Issues:**
- GitHub Issues: https://github.com/yourusername/EvoMotif/issues
- Include: Error message, code snippet, Python version
- Provide: Input data (if possible), expected behavior

---

## Appendix A: Mathematical Foundations

### Shannon Entropy Derivation

Shannon entropy quantifies uncertainty in a probability distribution:

$$H(X) = -\sum_{i=1}^{n} p(x_i) \log_2 p(x_i)$$

**Properties:**
- $H(X) = 0$ when one outcome has probability 1 (no uncertainty)
- $H(X)$ is maximized when all outcomes are equally likely
- Measured in bits (log base 2)

**Application to Conservation:**
- $X$ = amino acid at position $i$
- $p(x)$ = frequency of amino acid $x$
- $H = 0$ → perfectly conserved
- $H = \log_2(20) = 4.322$ → maximum variability

### BLOSUM62 Matrix

**How BLOSUM62 Was Created:**
1. Collected protein sequence blocks (aligned, conserved regions)
2. Counted observed amino acid pairs
3. Calculated expected frequencies under random model
4. Computed log-odds ratio:

$$S(a,b) = \lambda \log_2 \frac{P(a,b)}{P(a)P(b)}$$

Where:
- $P(a,b)$ = observed frequency of (a,b) pair
- $P(a)P(b)$ = expected frequency under independence
- $\lambda$ = scaling factor (2 bits)

**"62" in BLOSUM62:**
- Sequences >62% identical were clustered
- Representative sequence chosen from each cluster
- Removes bias from over-represented sequences

### Multiple Testing Correction

**Family-Wise Error Rate (FWER):**
$$P(\text{at least one false positive}) \leq \alpha$$

**False Discovery Rate (FDR):**
$$E\left[\frac{\text{False Positives}}{\text{Total Rejections}}\right] \leq \alpha$$

**Benjamini-Hochberg Procedure:**
1. Order p-values: $p_{(1)} \leq p_{(2)} \leq \cdots \leq p_{(m)}$
2. Find largest $k$ where: $p_{(k)} \leq \frac{k}{m} \alpha$
3. Reject null hypotheses for tests $1, 2, \ldots, k$

**Why FDR for Motif Discovery:**
- Expect multiple true positives (many conserved regions)
- More power than FWER control
- Controls expected proportion of false discoveries

---

## Appendix B: Glossary

**Conservation Score**: Quantitative measure of how similar amino acids are across species at a given position (0 = variable, 1 = conserved)

**Motif**: Short sequence pattern (typically 5-15 residues) that is conserved across multiple proteins

**BLOSUM62**: Substitution matrix capturing evolutionary relationships between amino acids, derived from aligned protein blocks

**Shannon Entropy**: Information-theoretic measure of variability/uncertainty in a distribution

**p-value**: Probability of observing data as extreme as the observed data, assuming the null hypothesis is true

**FDR**: False Discovery Rate - expected proportion of false positives among rejected null hypotheses

**pLDDT**: per-residue confidence score from AlphaFold, ranging 0-100 (>90 = high confidence)

**Effect Size**: Standardized measure of the magnitude of a difference, independent of sample size

**Gap Frequency**: Proportion of sequences with a gap (insertion/deletion) at a given alignment position

**Multiple Sequence Alignment (MSA)**: Alignment of three or more biological sequences (DNA, RNA, protein)

**Phylogenetic Tree**: Diagram showing evolutionary relationships among species or sequences

---

## Appendix C: References

**Key Publications:**

1. Henikoff S, Henikoff JG (1992). "Amino acid substitution matrices from protein blocks." PNAS 89(22):10915-9.

2. Shannon CE (1948). "A Mathematical Theory of Communication." Bell System Technical Journal 27:379-423.

3. Benjamini Y, Hochberg Y (1995). "Controlling the False Discovery Rate." Journal of the Royal Statistical Society B 57(1):289-300.

4. Katoh K, Standley DM (2013). "MAFFT Multiple Sequence Alignment Software." Molecular Biology and Evolution 30(4):772-780.

5. Price MN, Dehal PS, Arkin AP (2010). "FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments." PLoS ONE 5(3):e9490.

**Useful Resources:**

- NCBI Protein Database: https://www.ncbi.nlm.nih.gov/protein/
- PDB Structure Database: https://www.rcsb.org/
- AlphaFold Database: https://alphafold.ebi.ac.uk/
- MAFFT Documentation: https://mafft.cbrc.jp/alignment/software/
- Biopython Tutorial: https://biopython.org/DIST/docs/tutorial/Tutorial.html

---

**Document Version:** 1.0.0  
**Last Updated:** December 23, 2025  
**Authors:** EvoMotif Development Team  
**License:** MIT

---

*For questions, issues, or contributions, visit: https://github.com/yourusername/EvoMotif*
