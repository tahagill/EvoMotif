# EvoMotif: Evolutionary Protein Motif Discovery Tool

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-11%2F11-brightgreen.svg)](tests/)

**EvoMotif** is a comprehensive bioinformatics pipeline for discovering evolutionarily conserved protein motifs through multi-species sequence analysis, phylogenetic reconstruction, and 3D structural mapping.

---

## 🎯 Overview

EvoMotif integrates eight core modules to provide end-to-end analysis of protein conservation patterns:

1. **Sequence Retrieval** - Automated fetching from NCBI databases
2. **Multiple Sequence Alignment** - High-quality alignment with MAFFT
3. **Conservation Scoring** - Shannon entropy and BLOSUM62-based metrics
4. **Motif Discovery** - Sliding window and scattered residue detection
5. **Statistical Validation** - Permutation tests with FDR correction
6. **Phylogenetic Analysis** - Maximum likelihood tree inference
7. **3D Structure Mapping** - Conservation visualization on protein structures
8. **Results Compilation** - Comprehensive reports and visualizations

---

## 🚀 Quick Start

```bash
# Complete analysis with all 8 modules
python tests/run_complete_analysis.py AKT1 your.email@example.com \
    --pdb 3CQU \
    --max-sequences 50 \
    --min-conservation 0.70
```

**Output:** Interactive 3D structure, phylogenetic tree, conservation scores, and publication-ready figures!

---

## 📖 Documentation

| Document | Description |
|----------|-------------|
| **[PIPELINE_GUIDE.md](PIPELINE_GUIDE.md)** | Complete pipeline workflow and usage |
| **[STATISTICS_METHODS.md](STATISTICS_METHODS.md)** | Statistical methods with equations |
| **[PROJECT_REPORT.md](PROJECT_REPORT.md)** | Comprehensive project analysis |
| **[HOW_TO_VIEW_RESULTS.md](HOW_TO_VIEW_RESULTS.md)** | Visualization guide |

---

## 🔬 Validated Results

✅ **P53**: Found 7 conserved positions (5 Zn-binding cysteines, 1 DNA contact arginine)  
✅ **BRCA1**: Found 38 motifs (2 perfect tryptophans, 3 structural cysteines)  
✅ **AKT1**: Found DFG catalytic motif + activation loop (50 total motifs)

All results match known biological function! See [PROJECT_REPORT.md](PROJECT_REPORT.md) for details.

---

## 🛠️ Installation

```bash
# Python dependencies
pip install -r requirements.txt

# External tools (Ubuntu/Debian)
sudo apt-get install mafft fasttree
```

See full installation guide in [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md).

---

## 📊 Features

- ✅ Real data from NCBI (no fake/static data)
- ✅ Biologically accurate motif detection
- ✅ Statistical validation (permutation tests, FDR)
- ✅ Interactive 3D structure viewers
- ✅ Phylogenetic tree inference
- ✅ Publication-quality figures
- ✅ Fast execution (~30-60 sec per protein)

---

## 📄 Citation

```bibtex
@article{evomotif2025,
  title={EvoMotif: Evolutionary Protein Motif Discovery},
  author={Your Name},
  journal={Journal of Open Source Software},
  year={2025}
}
```

---

**Version:** 1.0.0 | **Status:** Production-ready | **Last Updated:** December 2, 2025
