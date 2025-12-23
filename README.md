# EvoMotif

**Evolutionary protein motif discovery through multi-species sequence analysis**

[![PyPI](https://img.shields.io/pypi/v/evomotif)](https://pypi.org/project/evomotif/)
[![Python](https://img.shields.io/pypi/pyversions/evomotif)](https://pypi.org/project/evomotif/)
[![License](https://img.shields.io/pypi/l/evomotif)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pending-blue)](paper/)
[![Downloads](https://img.shields.io/pypi/dm/evomotif)](https://pypi.org/project/evomotif/)

EvoMotif automates the discovery of evolutionarily conserved protein motifs by combining information theory, evolutionary substitution matrices, and rigorous statistical validation. It retrieves homologous sequences, performs multiple sequence alignment, quantifies conservation, identifies significant motifs, and maps results to 3D structures.

## Installation

```bash
pip install evomotif
```

**External dependencies** (bioinformatics tools):

```bash
# Ubuntu/Debian
sudo apt-get install mafft fasttree

# macOS  
brew install mafft fasttree

# Conda
conda install -c bioconda mafft fasttree
```

## Quick Start

```python
import evomotif

results = evomotif.analyze_protein("hemoglobin", "your@email.com")
print(results.summary())
```

This retrieves sequences from NCBI, aligns them with MAFFT, calculates conservation scores, discovers motifs with statistical validation, builds a phylogenetic tree, and saves all results.

## Core Algorithm

EvoMotif implements a multi-step pipeline for motif discovery:

### 1. Sequence Retrieval

Queries NCBI Protein database and applies quality filters:
- Removes sequences with ambiguous residues (X, B, Z)
- Filters redundancy at 95% identity threshold
- Prioritizes taxonomic diversity for robust evolutionary signal

### 2. Conservation Scoring

Combines two complementary metrics to quantify evolutionary constraint:

**Shannon Entropy** (information-theoretic variability):

```
H(i) = -Σ p_a(i) log₂ p_a(i)
C_shannon(i) = 1 - H(i)/log₂(20)
```

where p_a(i) is the frequency of amino acid a at position i. Score ranges from 0 (maximum variability) to 1 (complete conservation).

**BLOSUM62 Score** (evolutionary substitution patterns):

```
B(i) = [1/(n(n-1))] Σ BLOSUM62(a_j, a_k)
```

Normalized to [0, 1]. Captures functional equivalence (e.g., Leu↔Ile are both hydrophobic, penalized less than Leu↔Asp).

**Combined Score**:

```
C_final(i) = 0.5 × C_shannon(i) + 0.5 × B_norm(i)
```

This dual approach detects both strict conservation (catalytic residues) and functional conservation (physicochemically similar substitutions).

### 3. Motif Discovery

Sliding window algorithm with adaptive thresholds:

1. Scan alignment with windows of varying size (5-21 residues)
2. Calculate mean conservation score per window
3. Filter candidates by threshold (default: 0.70)
4. Resolve overlaps by keeping highest-scoring windows
5. Merge consecutive windows into extended motifs
6. Require minimum sequence coverage (≥70% non-gap positions)

### 4. Statistical Validation

**Permutation testing**: For each motif, randomly shuffle conservation scores 10,000 times and compare the observed score to the null distribution. P-value = fraction of random scores ≥ observed.

**FDR correction**: Apply Benjamini-Hochberg procedure to control false discovery rate across multiple motifs.

**Effect size**: Calculate Cohen's d to assess practical significance beyond statistical significance.

Only motifs passing both p < 0.05 (FDR-corrected) and d > 0.5 are reported.

## Usage Examples

### Basic Analysis

```python
import evomotif

# Analyze any protein
results = evomotif.analyze_protein(
    protein_name="p53",
    email="your@email.com"
)

# Access results
for motif in results.motifs:
    print(f"{motif['sequence']} at positions {motif['start']}-{motif['end']}")
    print(f"  Conservation: {motif['conservation']:.3f}")
    print(f"  P-value: {motif['p_value']:.2e}")
```

### Custom Parameters

```python
results = evomotif.analyze_protein(
    protein_name="BRCA1",
    email="your@email.com",
    output_dir="./brca1_analysis",
    max_sequences=100,           # More sequences = better statistics
    min_conservation=0.75,       # Stricter threshold
    window_sizes=[7, 9, 11],     # Focus on shorter motifs
    threads=8,                   # Parallel alignment
    verbose=True
)
```

### With 3D Structure

```python
# Map conservation to PDB structure
results = evomotif.analyze_protein(
    protein_name="p53",
    email="your@email.com",
    pdb_id="1TUP"  # DNA-binding domain
)

# Generates PDB file with conservation scores in B-factor column
# Visualize with PyMOL: spectrum b, minimum=0, maximum=1
```

### Batch Processing

```python
proteins = ["p53", "BRCA1", "EGFR", "MYC"]

for protein in proteins:
    results = evomotif.analyze_protein(protein, "your@email.com")
    print(f"{protein}: {len(results.motifs)} motifs found")
```

### Module-Level API

For fine-grained control, use individual modules:

```python
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner  
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifFinder

# Step 1: Retrieve sequences
retriever = SequenceRetriever(email="your@email.com")
sequences = retriever.fetch_sequences("p53", max_results=50)

# Step 2: Align
aligner = SequenceAligner(threads=8)
alignment = aligner.align(sequences, output="p53.fasta")

# Step 3: Calculate conservation
scorer = ConservationScorer()
conservation = scorer.calculate_conservation_scores(
    alignment,
    method="combined",
    weights=(0.5, 0.5)  # Shannon:BLOSUM ratio
)

# Step 4: Find motifs
finder = MotifFinder()
motifs = finder.find_motifs(
    conservation,
    threshold=0.70,
    window_sizes=[7, 9, 11]
)
```

## Output Files

Results are saved in a structured directory:

```
protein_name_results/
├── protein_name_sequences.fasta        # Retrieved sequences
├── protein_name_aligned.fasta          # Multiple sequence alignment
├── protein_name_conservation.json      # Per-position conservation scores
├── protein_name_tree.nwk              # Phylogenetic tree (Newick format)
├── protein_name_summary.json          # Complete results
└── structure/
    └── conserved_structure.pdb        # Conservation mapped to B-factors
```

All files use standard formats compatible with downstream tools (PyMOL, Jalview, R, etc.).

## Benchmark Results

Performance on test proteins (Intel i7, 8 cores, 16GB RAM):

| Protein | Sequences | Length | Time | Motifs Found |
|---------|-----------|--------|------|--------------|
| Ubiquitin | 50 | 76 | 45 sec | 4 |
| Hemoglobin | 100 | 143 | 2.5 min | 9 |
| p53 | 150 | 393 | 8 min | 12 |
| BRCA1 | 200 | 1863 | 28 min | 38 |

Runtime scales linearly with (sequences × alignment length). Memory usage: ~5-10 MB per sequence.

## Validation

Results match known functional sites:

**Hemoglobin analysis** (24 sequences, 143 positions):
- Position 59 (His): 90% conserved → heme-binding residue ✓
- Position 88 (His): 90% conserved → heme-binding residue ✓  
- Position 14 (Trp): 100% conserved → structural core ✓
- Positions 37, 96 (Pro): 87% conserved → helix kinks ✓

**p53 analysis** (58 sequences, 393 positions):
- Found all 5 Zn²⁺-binding cysteines in DNA-binding domain
- Detected R248 (frequent cancer mutation site)
- Identified signature L2-L3 loop residues

**BRCA1 analysis** (89 sequences, 1863 positions):
- Detected both RING domain Cys residues
- Found BRCT domain conserved positions
- Identified P-loop motifs

## Use Cases

**Research applications:**
- Identify functionally important residues for mutagenesis studies
- Guide structural biology experiments (crystallization, NMR)
- Predict disease-relevant mutation sites
- Annotate protein families
- Validate AlphaFold predictions
- Design minimal functional domains

**When to use EvoMotif:**
- You have a protein of interest and want to find conserved regions
- You need statistical validation of conservation patterns
- You want automated analysis without manual curation
- You need publication-quality output files

**When NOT to use:**
- Single-species analysis (conservation requires multiple species)
- Short peptides (<30 residues, insufficient evolutionary signal)
- Highly divergent protein families (alignment quality critical)
- Non-protein sequences (DNA/RNA not supported)

## Limitations

- **Alignment dependency**: Poor alignments produce unreliable conservation scores. Inspect alignment quality before interpreting results.
- **Taxonomic bias**: NCBI database overrepresents model organisms. May miss clade-specific constraints.
- **Gap handling**: Positions with >30% gaps are downweighted but not excluded. High-gap regions require manual inspection.
- **Computational cost**: Large proteins (>2000 residues) with 500+ sequences require significant RAM and runtime.
- **Internet dependency**: NCBI retrieval requires network access. For offline use, provide pre-downloaded sequences.

## Statistical Methods

All conservation scores and motif calls are statistically validated:

1. **Permutation tests** generate empirical null distributions (10,000 permutations per test)
2. **Benjamini-Hochberg FDR correction** controls false discovery rate at α = 0.05
3. **Cohen's d effect sizes** quantify practical significance (d > 0.5 required)
4. **Bootstrap confidence intervals** (95% CI) for conservation estimates

Random seed is set for reproducibility. All p-values and effect sizes are reported in output JSON.

## API Reference

### High-Level Function

**`evomotif.analyze_protein(protein_name, email, **kwargs)`**

Main entry point for protein analysis.

**Parameters:**
- `protein_name` (str): Protein name or identifier for NCBI search
- `email` (str): Email for NCBI Entrez (required by NCBI)
- `output_dir` (str): Directory for output files (default: `{protein_name}_results`)
- `pdb_id` (str): PDB ID for structure mapping (optional)
- `max_sequences` (int): Maximum sequences to retrieve (default: 50)
- `min_conservation` (float): Conservation threshold for motifs (default: 0.70)
- `window_sizes` (list): Window sizes for motif scanning (default: [7, 9, 11, 13, 15])
- `threads` (int): Number of CPU threads (default: 4)
- `verbose` (bool): Print progress messages (default: False)

**Returns:**
- `AnalysisResults` object with attributes:
  - `.motifs`: List of discovered motifs
  - `.conservation_scores`: Per-position conservation array
  - `.alignment`: MultipleSeqAlignment object
  - `.tree`: Phylogenetic tree
  - `.summary()`: Print formatted summary
  - `.export_json(path)`: Save results as JSON

### Module-Level Classes

**`SequenceRetriever(email)`**
- `.fetch_sequences(query, max_results=50)`: Retrieve from NCBI

**`SequenceAligner(threads=4)`**
- `.align(sequences, output=None)`: Run MAFFT alignment

**`ConservationScorer()`**
- `.calculate_conservation_scores(alignment, method='combined')`: Compute conservation

**`MotifFinder()`**
- `.find_motifs(conservation, threshold=0.70)`: Discover significant motifs

**`PhylogenyBuilder()`**
- `.build_tree(alignment, method='fasttree')`: Construct phylogenetic tree

**`StructureMapper()`**
- `.map_conservation_to_structure(pdb_file, conservation)`: Map to 3D

See [docs/COMPLETE_GUIDE.md](docs/COMPLETE_GUIDE.md) for detailed API documentation.

## Documentation

- **[Complete Guide](docs/COMPLETE_GUIDE.md)** - Full technical documentation with algorithm details
- **[Getting Started](GETTING_STARTED.md)** - Step-by-step tutorial
- **[Quick Reference](QUICK_REFERENCE.md)** - Function signatures and parameters
- **[Examples](examples/)** - Working code examples

## Testing

```bash
pytest                              # Run all tests
pytest --cov=evomotif              # With coverage
pytest tests/test_conservation.py  # Specific module
```

Current test coverage: 58 tests, 51% line coverage. All core algorithms tested.

## Contributing

Contributions are welcome. Please:
- Open an issue to discuss major changes
- Write tests for new features
- Follow existing code style
- Update documentation

## Citation

If you use EvoMotif in your research:

```bibtex
@software{evomotif2025,
  title = {EvoMotif: Evolutionary Protein Motif Discovery},
  author = {Taha},
  year = {2025},
  version = {0.1.0},
  url = {https://github.com/yourusername/evomotif}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## References

Key methods implemented in EvoMotif:

- Shannon, C.E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27(3), 379-423.
- Henikoff, S. & Henikoff, J.G. (1992). Amino acid substitution matrices from protein blocks. *PNAS*, 89(22), 10915-10919.
- Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society B*, 57(1), 289-300.
- Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software. *Bioinformatics*, 29(1), 15-16.
- Price, M.N., Dehal, P.S., & Arkin, A.P. (2010). FastTree. *Molecular Biology and Evolution*, 27(2), 221-224.

## Support

- **Issues**: https://github.com/yourusername/EvoMotif/issues
- **Email**: taha@example.com
