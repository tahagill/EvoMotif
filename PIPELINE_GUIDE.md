# EvoMotif Pipeline Guide

Complete workflow documentation for the EvoMotif protein motif discovery pipeline.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Installation](#installation)
3. [Module Descriptions](#module-descriptions)
4. [Usage Examples](#usage-examples)
5. [Input/Output Specifications](#inputoutput-specifications)
6. [Configuration Options](#configuration-options)
7. [Troubleshooting](#troubleshooting)

---

## Pipeline Overview

### Complete Workflow

```
[1] Sequence Retrieval (NCBI)
         ↓
[2] Multiple Alignment (MAFFT)
         ↓
[3] Conservation Scoring
         ↓
[4] Motif Discovery
         ↓
[5] Statistical Validation
         ↓
[6] Phylogenetic Analysis
         ↓
[7] 3D Structure Mapping
         ↓
[8] Results Compilation
```

### Processing Time

| Step | Typical Duration | Notes |
|------|-----------------|-------|
| Sequence Retrieval | 2-5 seconds | Depends on NCBI server |
| Alignment | 5-30 seconds | Scales with sequence count |
| Conservation | 1-3 seconds | Fast computation |
| Motif Discovery | 1-2 seconds | Both methods combined |
| Statistics | 2-10 seconds | Permutation tests |
| Phylogeny | 2-5 seconds | FastTree is fast! |
| Structure Mapping | 1-3 seconds | PDB download time |
| Total | 15-60 seconds | End-to-end |

---

## Installation

### Prerequisites

**Operating System:**
- Linux (Ubuntu 20.04+, Debian, CentOS)
- macOS (10.14+)
- Windows (WSL2 recommended)

**Python Version:**
- Python 3.8 or higher

### Step-by-Step Installation

#### 1. Clone Repository

```bash
git clone https://github.com/yourusername/evomotif.git
cd evomotif
```

#### 2. Create Virtual Environment

```bash
python3 -m venv .venv
source .venv/bin/activate  # Linux/Mac
# or
.venv\Scripts\activate  # Windows
```

#### 3. Install Python Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

**Required packages:**
```
biopython>=1.79
numpy>=1.21.0
scipy>=1.7.0
pandas>=1.3.0
matplotlib>=3.4.0
seaborn>=0.11.0
scikit-learn>=0.24.0
statsmodels>=0.12.0
ete3>=3.1.2
py3Dmol>=2.0.0
requests>=2.26.0
```

#### 4. Install External Tools

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y mafft fasttree
```

**macOS:**
```bash
brew install mafft fasttree
```

**Manual Installation:**
- MAFFT: https://mafft.cbrc.jp/alignment/software/
- FastTree: http://www.microbesonline.org/fasttree/

#### 5. Verify Installation

```bash
# Check Python packages
python -c "import Bio; print(f'BioPython: {Bio.__version__}')"

# Check external tools
mafft --version
fasttree 2>&1 | head -1

# Run tests
python -m pytest tests/test_integration.py -v
```

---

## Module Descriptions

### Module 1: Sequence Retrieval

**Purpose:** Fetch protein sequences from NCBI protein database

**Implementation:** `evomotif/retrieval.py`

**Key Functions:**
- `retrieve_all()` - Complete retrieval pipeline
- `search_protein()` - Search NCBI with query
- `fetch_sequences()` - Download sequences in batches
- `filter_sequences()` - Quality control filtering
- `remove_duplicates()` - Deduplicate by sequence identity

**Usage:**
```python
from evomotif.retrieval import SequenceRetriever

retriever = SequenceRetriever(email="your@email.com")
sequences = retriever.retrieve_all(
    protein_name="AKT1",
    output_path="sequences.fasta",
    max_sequences=50,
    min_length=100,
    max_length=3000,
    remove_fragments=False
)
```

**Parameters:**
- `protein_name`: Protein name for NCBI search
- `output_path`: Output FASTA file path
- `max_sequences`: Maximum sequences to retrieve (default: 1000)
- `min_length`: Minimum sequence length (default: 30)
- `max_length`: Maximum sequence length (default: None)
- `remove_fragments`: Remove partial sequences (default: False)

**Output:** FASTA file with retrieved sequences

---

### Module 2: Multiple Sequence Alignment

**Purpose:** Create high-quality multiple sequence alignment

**Implementation:** `evomotif/alignment.py`

**Key Functions:**
- `align_sequences()` - MAFFT alignment wrapper
- `calculate_alignment_stats()` - Alignment quality metrics

**Usage:**
```python
from evomotif.alignment import SequenceAligner

aligner = SequenceAligner()
alignment = aligner.align_sequences(
    input_fasta="sequences.fasta",
    output_fasta="aligned.fasta",
    method="linsi",  # L-INS-i algorithm
    threads=4
)

stats = aligner.calculate_alignment_stats(alignment)
```

**MAFFT Methods:**
- `linsi`: L-INS-i (accurate, slow) - **DEFAULT**
- `ginsi`: G-INS-i (global alignment)
- `einsi`: E-INS-i (conserved domains)
- `auto`: Automatic selection

**Alignment Statistics:**
- `n_sequences`: Number of sequences
- `alignment_length`: Total positions
- `gap_percentage`: Overall gap frequency
- `gappy_columns`: Columns with >50% gaps
- `gappy_columns_pct`: Percentage of gappy columns

**Output:** Aligned FASTA file

---

### Module 3: Conservation Scoring

**Purpose:** Calculate position-wise conservation scores

**Implementation:** `evomotif/conservation.py`

**Key Functions:**
- `calculate_combined_conservation()` - Combined metric
- `calculate_shannon_entropy()` - Information content
- `calculate_blosum_score()` - Substitution scoring
- `calculate_gap_frequency()` - Gap distribution
- `get_consensus_sequence()` - Consensus extraction

**Usage:**
```python
from evomotif.conservation import ConservationScorer

scorer = ConservationScorer()

# Combined conservation (Shannon + BLOSUM62)
conservation = scorer.calculate_combined_conservation(alignment)

# Individual metrics
shannon = scorer.calculate_shannon_entropy(alignment)
blosum = scorer.calculate_blosum_score(alignment)
gaps = scorer.calculate_gap_frequency(alignment)
consensus = scorer.get_consensus_sequence(alignment)
```

**Conservation Metrics:**

1. **Shannon Entropy** (H):
   - Range: 0 (conserved) to high (variable)
   - Normalized to 0-1 scale

2. **BLOSUM62 Score**:
   - Substitution matrix-based
   - Higher = more conserved

3. **Combined Score**:
   - Weighted average of Shannon + BLOSUM
   - Range: 0 (variable) to 1 (conserved)

**Output:** Conservation scores JSON file

---

### Module 4: Motif Discovery

**Purpose:** Identify conserved motifs using two complementary methods

**Implementation:** `evomotif/motif_discovery.py`

**Key Functions:**
- `find_motifs()` - Sliding window approach (consecutive motifs)
- `find_scattered_motifs()` - Position-wise approach (scattered residues)
- `merge_overlapping_motifs()` - Combine adjacent motifs
- `export_motif_sequences()` - Extract motif FASTA files

**Usage:**
```python
from evomotif.motif_discovery import MotifDiscoverer

discoverer = MotifDiscoverer()

# Method 1: Consecutive motifs (sliding window)
motifs = discoverer.find_motifs(
    alignment,
    conservation_scores,
    min_conservation=0.7,
    min_length=5,
    max_gap_frequency=0.3
)

# Method 2: Scattered conserved residues
scattered = discoverer.find_scattered_motifs(
    conservation_scores,
    gap_frequencies,
    consensus_sequence,
    min_conservation=0.7,
    max_gap=0.5
)
```

**Motif Discovery Algorithms:**

**Algorithm 1: Sliding Window**
1. Scan alignment with window (default: 5 positions)
2. Calculate average conservation in window
3. Filter by conservation threshold
4. Merge overlapping windows
5. Output consecutive motifs

**Algorithm 2: Scattered Positions**
1. Evaluate each position independently
2. Filter by conservation + gap threshold
3. Identify individual conserved residues
4. No requirement for adjacency

**When to Use Each:**
- **Sliding window**: Domain motifs, linear epitopes, binding sites
- **Scattered**: Metal-binding sites, catalytic triads, structural residues

**Output:** Motif metadata JSON + individual FASTA files

---

### Module 5: Statistical Validation

**Purpose:** Validate motif significance with permutation tests

**Implementation:** `evomotif/stats.py`

**Key Functions:**
- `permutation_test()` - Non-parametric significance testing
- `multiple_testing_correction()` - FDR correction (Benjamini-Hochberg)
- `calculate_effect_size()` - Cohen's d
- `bootstrap_confidence_interval()` - Bootstrap CI estimation

**Usage:**
```python
from evomotif.stats import StatisticalValidator

validator = StatisticalValidator()

# Permutation test for each motif
for motif in motifs:
    result = validator.permutation_test(
        observed_score=motif.conservation,
        null_distribution=background_conservation,
        n_permutations=1000,
        alternative='greater'
    )
    
# FDR correction
p_values = [motif.p_value for motif in motifs]
corrected = validator.multiple_testing_correction(
    p_values,
    method='fdr_bh',  # Benjamini-Hochberg
    alpha=0.05
)
```

**Statistical Methods:**
- Permutation tests (non-parametric)
- Benjamini-Hochberg FDR correction
- Cohen's d effect size
- Bootstrap confidence intervals

See [STATISTICS_METHODS.md](STATISTICS_METHODS.md) for detailed equations.

**Output:** Updated motif metadata with p-values and FDR

---

### Module 6: Phylogenetic Analysis

**Purpose:** Infer evolutionary relationships among sequences

**Implementation:** `evomotif/phylogeny.py`

**Key Functions:**
- `build_tree_fasttree()` - Maximum likelihood tree inference
- `save_tree()` - Export Newick format
- `calculate_tree_stats()` - Tree metrics

**Usage:**
```python
from evomotif.phylogeny import PhylogeneticAnalyzer

analyzer = PhylogeneticAnalyzer()

tree = analyzer.build_tree_fasttree(
    alignment_file="aligned.fasta",
    output_file="tree.nwk",
    model="JTT"  # Amino acid substitution model
)
```

**Substitution Models:**
- `JTT`: Jones-Taylor-Thornton (general purpose) - **DEFAULT**
- `WAG`: Whelan and Goldman (broad coverage)
- `LG`: Le and Gascuel (recent model)

**Tree Output:**
- Newick format (.nwk)
- Branch lengths in substitutions/site
- Bootstrap support (if requested)

**Output:** Phylogenetic tree in Newick format

---

### Module 7: 3D Structure Mapping

**Purpose:** Map conservation onto protein 3D structures

**Implementation:** `evomotif/structure.py`

**Key Functions:**
- `load_structure()` - Load PDB file
- `map_alignment_to_structure()` - Position mapping
- `visualize_structure_3d()` - Interactive py3Dmol viewer

**Usage:**
```python
from evomotif.structure import StructureMapper

mapper = StructureMapper()

# Load PDB structure
structure = mapper.load_structure("3CQU.pdb")

# Map alignment positions to structure
mapping = mapper.map_alignment_to_structure(
    structure,
    alignment_sequence,
    chain_id='A'
)

# Create 3D visualization
html = mapper.visualize_structure_3d(
    pdb_file="3CQU.pdb",
    motif_residues=[234, 250, 345, 631, 710],
    output_html="structure_3d.html"
)
```

**Visualization Features:**
- Interactive rotation/zoom
- Conserved residues highlighted (RED sticks)
- Protein backbone (GRAY cartoon)
- Self-contained HTML (no internet required)

**Output:** Interactive HTML file with embedded 3D viewer

---

### Module 8: Results Compilation

**Purpose:** Generate comprehensive analysis reports

**Implementation:** Integrated in `run_complete_analysis.py`

**Output Files:**

```
complete_analysis_results/{PROTEIN}/
├── {PROTEIN}_sequences.fasta       # Retrieved sequences
├── {PROTEIN}_aligned.fasta         # Multiple alignment
├── {PROTEIN}_conservation.json     # All conservation scores
├── conserved_positions.json        # Scattered conserved residues
├── {PROTEIN}_tree.nwk             # Phylogenetic tree
├── {PROTEIN}_structure_3d.html    # Interactive 3D viewer
├── {PDB_ID}.pdb                   # PDB structure file
├── {PROTEIN}_summary.json         # Complete summary
└── consecutive_motifs/            # Individual motif FASTA files
    ├── motif_1.fasta
    ├── motif_2.fasta
    └── ...
```

---

## Usage Examples

### Example 1: Basic Analysis

```bash
python tests/run_complete_analysis.py p53 your@email.com \
    --max-sequences 30 \
    --min-conservation 0.70
```

### Example 2: With Structure Visualization

```bash
python tests/run_complete_analysis.py AKT1 your@email.com \
    --pdb 3CQU \
    --max-sequences 50 \
    --min-conservation 0.70
```

### Example 3: Custom Parameters

```bash
python tests/run_complete_analysis.py BRCA1 your@email.com \
    --pdb 1JM7 \
    --max-sequences 100 \
    --min-conservation 0.65 \
    --threads 8 \
    --output custom_results/ \
    --verbose
```

### Example 4: Module-by-Module

```python
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer
from pathlib import Path

# Step 1: Retrieve
retriever = SequenceRetriever(email="your@email.com")
sequences = retriever.retrieve_all(
    protein_name="ubiquitin",
    output_path=Path("ubiquitin_seqs.fasta"),
    max_sequences=20
)

# Step 2: Align
aligner = SequenceAligner()
alignment = aligner.align_sequences(
    input_fasta="ubiquitin_seqs.fasta",
    output_fasta="ubiquitin_aligned.fasta"
)

# Step 3: Score conservation
scorer = ConservationScorer()
conservation = scorer.calculate_combined_conservation(alignment)

# Step 4: Find motifs
discoverer = MotifDiscoverer()
motifs = discoverer.find_motifs(
    alignment,
    conservation,
    min_conservation=0.75
)

print(f"Found {len(motifs)} conserved motifs!")
```

---

## Input/Output Specifications

### Input Requirements

**Protein Name:**
- Standard gene/protein name (case-insensitive)
- Examples: "p53", "BRCA1", "AKT1", "ubiquitin"

**Email (Required for NCBI):**
- Valid email address
- Used for NCBI Entrez API identification
- Not shared or stored

**PDB ID (Optional):**
- 4-character PDB identifier
- Example: "1TUP", "3CQU", "1JM7"
- If not provided, structure mapping skipped

### Output File Formats

**FASTA (.fasta)**
```
>NP_000537.3 tumor protein p53 [Homo sapiens]
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPG
...
```

**JSON (.json)**
```json
{
  "position": 215,
  "conservation": 0.877,
  "gap": 0.379,
  "amino_acid": "C"
}
```

**Newick (.nwk)**
```
((seq1:0.1,seq2:0.15):0.2,(seq3:0.05,seq4:0.12):0.18);
```

**HTML (.html)**
- Self-contained interactive 3D viewer
- Embedded py3Dmol library
- Works offline

---

## Configuration Options

### Command-Line Arguments

```
python run_complete_analysis.py PROTEIN EMAIL [OPTIONS]

Positional Arguments:
  PROTEIN              Protein name to analyze
  EMAIL                Email for NCBI Entrez API

Optional Arguments:
  --pdb PDB_ID         PDB structure ID for 3D visualization
  --max-sequences N    Maximum sequences to retrieve (default: 50)
  --min-conservation F Conservation threshold 0-1 (default: 0.70)
  --threads N          CPU threads for alignment (default: 4)
  --output DIR         Output directory (default: complete_analysis_results)
  --verbose            Enable debug logging
  --help               Show help message
```

### Python API Configuration

```python
# Conservation thresholds
CONSERVATION_THRESHOLD = 0.70
MIN_MOTIF_LENGTH = 5
MAX_GAP_FREQUENCY = 0.30

# Statistical parameters
N_PERMUTATIONS = 1000
FDR_ALPHA = 0.05

# Alignment parameters
MAFFT_METHOD = "linsi"
ALIGNMENT_THREADS = 4

# Phylogeny parameters
TREE_MODEL = "JTT"
BOOTSTRAP_REPLICATES = 0  # Set to 100 for bootstrap support
```

---

## Troubleshooting

### Common Issues

**1. "MAFFT not found"**
```bash
# Install MAFFT
sudo apt-get install mafft

# Or check PATH
which mafft
```

**2. "FastTree not found"**
```bash
# Install FastTree
sudo apt-get install fasttree

# Or check PATH
which FastTree
```

**3. "No sequences found"**
- Check protein name spelling
- Try broader search terms
- Check NCBI server status
- Verify email is valid

**4. "Too few sequences"**
- Increase `--max-sequences`
- Adjust length filters (min_length, max_length)
- Set `remove_fragments=False`

**5. "Structure mapping failed"**
- Verify PDB ID is correct
- Check internet connection (for PDB download)
- Ensure sequence matches structure

**6. "Import errors"**
```bash
# Reinstall dependencies
pip install --force-reinstall -r requirements.txt
```

### Performance Optimization

**Large datasets (>100 sequences):**
```bash
# Use more threads
--threads 8

# Use faster MAFFT method
# Edit: alignment.py, change method to "auto"

# Skip bootstrap (phylogeny)
# Faster tree inference
```

**Memory issues:**
```bash
# Reduce max sequences
--max-sequences 50

# Use less RAM-intensive methods
```

### Logging and Debugging

```bash
# Enable verbose output
python run_complete_analysis.py PROTEIN EMAIL --verbose

# Check log files
tail -f evomotif.log

# Python debugging
python -m pdb run_complete_analysis.py PROTEIN EMAIL
```

---

## Best Practices

### 1. Email Configuration
- Use institutional email for NCBI access
- Set environment variable: `export NCBI_EMAIL=your@email.com`

### 2. Sequence Retrieval
- Start with 20-50 sequences
- Increase if needed for statistical power
- Balance: More sequences = better statistics, longer runtime

### 3. Conservation Thresholds
- **Conservative** (0.85+): High-confidence motifs only
- **Moderate** (0.70-0.85): Balance sensitivity/specificity
- **Permissive** (0.60-0.70): Exploratory analysis

### 4. Statistical Validation
- Always use FDR correction for multiple testing
- Report both raw and corrected p-values
- Consider effect sizes, not just p-values

### 5. Visualization
- Use interactive HTML for exploration
- Generate PNG/PDF for publications
- Upload trees to iTOL for beautiful figures

### 6. Results Interpretation
- Cross-validate with literature
- Check protein databases (UniProt, PDB)
- Verify with experimental data when available

---

## Next Steps

- See [STATISTICS_METHODS.md](STATISTICS_METHODS.md) for statistical details
- See [PROJECT_REPORT.md](PROJECT_REPORT.md) for complete analysis examples
- See [HOW_TO_VIEW_RESULTS.md](HOW_TO_VIEW_RESULTS.md) for visualization guide

---

**Last Updated:** December 2, 2025  
**Version:** 1.0.0
