# Getting Started with EvoMotif

Welcome! This guide will get you analyzing proteins in **under 5 minutes**.

---

## 🚀 Quick Installation

```bash
# Clone repository
git clone https://github.com/tahagill/EvoMotif.git
cd EvoMotif

# Install dependencies
pip install -r requirements.txt

# Install external tools (Ubuntu/Debian)
sudo apt-get install mafft fasttree

# Or on macOS
brew install mafft fasttree

# Or with Conda
conda install -c bioconda mafft fasttree
```

---

## 🎯 Your First Analysis (2 Lines!)

```python
import evomotif

results = evomotif.analyze_protein("p53", "your@email.com")
print(results.summary())
```

**Output:**
```
==============================================================
EvoMotif Analysis Results: p53
==============================================================
Sequences analyzed: 50
Consecutive motifs found: 12
Conserved positions: 87
Mean conservation: 0.654
Max conservation: 0.982

Top 5 motifs:
  1. CQKFSY (pos 176-181, cons=0.892)
  2. TFRHSV (pos 245-250, cons=0.865)
  ...

📁 Results saved to: ./evomotif_results/p53
==============================================================
```

**That's it!** 🎉 You just discovered evolutionarily conserved motifs in p53.

---

## 📊 What Just Happened?

Behind the scenes, EvoMotif:

1. ✅ **Retrieved 50 p53 sequences** from NCBI across multiple species
2. ✅ **Aligned them** using MAFFT L-INS-i algorithm
3. ✅ **Calculated conservation** using Shannon entropy + BLOSUM62
4. ✅ **Discovered motifs** using sliding windows + scattered residues
5. ✅ **Validated statistically** with permutation tests & FDR correction
6. ✅ **Built phylogenetic tree** using FastTree maximum likelihood
7. ✅ **Saved everything** to `evomotif_results/p53/`

---

## 🧬 Example 2: With 3D Structure

```python
import evomotif

# Add PDB structure visualization
results = evomotif.analyze_protein(
    protein_name="p53",
    email="your@email.com",
    pdb_id="1TUP"  # P53 DNA-binding domain structure
)

print(results.summary())

# Access 3D structure file
pdb_file = results.get_file('pdb')
print(f"Structure saved to: {pdb_file}")
```

---

## ⚙️ Example 3: Custom Parameters

```python
import evomotif

results = evomotif.analyze_protein(
    protein_name="BRCA1",
    email="your@email.com",
    output_dir="./my_results/brca1",  # Custom output location
    max_sequences=100,                 # More sequences = better signal
    min_conservation=0.65,             # Lower threshold = more motifs
    threads=8,                         # Use more CPU cores
    verbose=True                       # Show detailed logging
)

# Access specific results
print(f"Found {len(results.motifs)} consecutive motifs")
print(f"Found {len(results.conserved_positions)} conserved positions")

# Get top 10 conserved positions
sorted_pos = sorted(
    results.conserved_positions,
    key=lambda x: x['conservation'],
    reverse=True
)

for i, pos in enumerate(sorted_pos[:10], 1):
    print(f"{i:2d}. Position {pos['position']:3d}: {pos['amino_acid']} "
          f"(conservation={pos['conservation']:.3f})")
```

---

## 📁 Understanding the Output

After running analysis, you'll find these files in the output directory:

```
evomotif_results/p53/
├── p53_sequences.fasta              # Retrieved sequences
├── p53_aligned.fasta                # Multiple sequence alignment
├── p53_conservation.json            # Conservation scores for all positions
├── p53_tree.nwk                     # Phylogenetic tree (Newick format)
├── p53_summary.json                 # Complete results summary
├── conserved_positions.json         # Scattered conserved residues
└── motifs/                          # Individual motif FASTA files
    ├── motif_1.fasta
    ├── motif_2.fasta
    └── ...
```

---

## 🔍 Checking Dependencies

Before running large analyses, check that tools are installed:

```python
import evomotif

pipeline = evomotif.EvoMotifPipeline()
deps = pipeline.check_dependencies()

for tool, available in deps.items():
    status = "✓" if available else "✗"
    print(f"{status} {tool}")
```

**Output:**
```
✓ mafft
✓ fasttree
```

If you see ✗, install the missing tool:
```bash
# Ubuntu/Debian
sudo apt-get install mafft fasttree

# macOS
brew install mafft fasttree

# Conda
conda install -c bioconda mafft fasttree
```

---

## 💾 Saving and Loading Results

```python
import evomotif
from pathlib import Path

# Run analysis
results = evomotif.analyze_protein("ubiquitin", "your@email.com")

# Export to JSON
results.export_json(Path("my_results.json"))

# Access files later
alignment_file = results.get_file('alignment')
conservation_file = results.get_file('conservation')
tree_file = results.get_file('tree')
```

---

## 🎓 Next Steps

### For Scientists (Use Simple API)

Just use `evomotif.analyze_protein()` for everything! See:
- **[examples/simple_api_demo.py](examples/simple_api_demo.py)** - 6 complete examples
- **[USER_GUIDE.md](USER_GUIDE.md)** - Comprehensive guide with theory

### For Developers (Use Individual Modules)

Want fine-grained control? Use individual modules:
- **[PIPELINE_GUIDE.md](PIPELINE_GUIDE.md)** - Module-by-module documentation
- **[examples/01_simple_analysis.py](examples/01_simple_analysis.py)** - Step-by-step example

### Advanced Topics

- **[STATISTICS_METHODS.md](STATISTICS_METHODS.md)** - Statistical validation details
- **[ALPHAFOLD_INTEGRATION.md](ALPHAFOLD_INTEGRATION.md)** - AlphaFold confidence analysis
- **[HOW_TO_VIEW_RESULTS.md](HOW_TO_VIEW_RESULTS.md)** - Visualization guide

---

## 🐛 Troubleshooting

### "MAFFT not found"
```bash
sudo apt-get install mafft  # Ubuntu
brew install mafft          # macOS
```

### "FastTree not found"
```bash
sudo apt-get install fasttree  # Ubuntu
brew install fasttree          # macOS
```

### "Too few sequences"
Try increasing `max_sequences`:
```python
results = evomotif.analyze_protein(
    "your_protein",
    "your@email.com",
    max_sequences=100  # Increase from default 50
)
```

### "No motifs found"
Try lowering conservation threshold:
```python
results = evomotif.analyze_protein(
    "your_protein",
    "your@email.com",
    min_conservation=0.60  # Lower from default 0.70
)
```

### Analysis is slow
Use fewer sequences and more threads:
```python
results = evomotif.analyze_protein(
    "your_protein",
    "your@email.com",
    max_sequences=30,  # Fewer sequences
    threads=8          # More CPU cores
)
```

---

## 📚 API Reference

### `analyze_protein()`

```python
evomotif.analyze_protein(
    protein_name: str,          # Required: protein name (e.g., "p53")
    email: str,                 # Required: your email for NCBI
    output_dir: str = None,     # Optional: output directory
    pdb_id: str = None,         # Optional: PDB ID for structure
    max_sequences: int = 50,    # Optional: max sequences to retrieve
    min_conservation: float = 0.70,  # Optional: conservation threshold
    threads: int = 4,           # Optional: CPU threads
    verbose: bool = False,      # Optional: debug logging
    check_deps: bool = True     # Optional: check dependencies
) -> AnalysisResults
```

### `AnalysisResults` Object

```python
results.protein                  # Protein name
results.n_sequences             # Number of sequences
results.motifs                  # List of consecutive motifs
results.conserved_positions     # List of conserved positions
results.conservation_scores     # Conservation array
results.output_dir              # Output directory path
results.files                   # Dict of generated files

results.summary()               # Print human-readable summary
results.get_file('alignment')   # Get specific file path
results.export_json(path)       # Export to JSON
```

---

## ✨ Pro Tips

1. **Start with well-known proteins** (p53, ubiquitin, insulin) to validate
2. **Use 30-50 sequences** for good balance of speed vs. accuracy
3. **Lower conservation threshold** (0.60-0.65) for divergent protein families
4. **Add PDB structure** when available for 3D visualization
5. **Check conserved positions** even if no consecutive motifs found
6. **Use verbose=True** if analysis fails to see detailed errors

---

## 🎉 Success!

You're now ready to discover conserved motifs in any protein! 

**Questions?** See [USER_GUIDE.md](USER_GUIDE.md) for detailed documentation.

**Issues?** Open an issue on [GitHub](https://github.com/tahagill/EvoMotif/issues).

Happy analyzing! 🧬
