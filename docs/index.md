# EvoMotif Documentation

Welcome to the EvoMotif documentation!

## Overview

EvoMotif is an evolution-driven framework for discovering novel protein motifs, mapping them to 3D structures, and analyzing their evolutionary origin.

## Features

- **Automated sequence retrieval** from NCBI databases
- **Multiple sequence alignment** using MAFFT
- **Conservation scoring** with Shannon entropy and BLOSUM62
- **Motif discovery** using sliding-window analysis
- **HMM profiling** for motif generalization
- **Phylogenetic analysis** to determine motif evolutionary origin
- **3D structural mapping** with py3Dmol
- **Variant analysis** for pathogenic mutation enrichment
- **Comprehensive statistical testing**

## Quick Start

```python
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer

# Retrieve sequences
retriever = SequenceRetriever(email="your@email.com")
sequences = retriever.retrieve_all(
    "p53",
    output_path="data/p53.fasta",
    max_sequences=500
)

# Align sequences
aligner = SequenceAligner()
alignment = aligner.align_sequences(
    input_fasta="data/p53.fasta",
    output_fasta="data/p53_aligned.fasta",
    method="linsi"
)

# Calculate conservation
scorer = ConservationScorer()
conservation = scorer.calculate_combined_conservation(alignment)

# Discover motifs
discoverer = MotifDiscoverer()
motifs = discoverer.discover_motifs(
    alignment,
    conservation,
    gap_frequency,
    consensus
)
```

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/EvoMotif.git
cd EvoMotif

# Install dependencies
pip install -r requirements.txt

# Install EvoMotif
pip install -e .
```

### External Dependencies

EvoMotif requires the following external tools:

- **MAFFT** (v7.0+): Multiple sequence alignment
- **HMMER** (v3.0+): HMM profiling
- **FastTree** or **IQ-TREE**: Phylogenetic tree construction
- **DSSP** (optional): Secondary structure assignment

See [Installation Guide](installation.md) for detailed instructions.

## Table of Contents

```{toctree}
:maxdepth: 2

installation
tutorials/index
api/index
examples/index
contributing
```

## Citation

If you use EvoMotif in your research, please cite:

```bibtex
@article{evomotif2025,
  title={EvoMotif: An Evolution-Driven Framework for Discovering Novel Protein Motifs},
  author={Taha},
  journal={Journal of Open Source Software},
  year={2025}
}
```

## License

EvoMotif is licensed under the MIT License. See [LICENSE](../LICENSE) for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/EvoMotif/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/EvoMotif/discussions)
- **Email**: your@email.com
