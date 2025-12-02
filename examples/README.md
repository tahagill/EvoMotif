# Example Scripts

This directory contains example scripts demonstrating EvoMotif functionality.

## Basic Examples

- `01_simple_analysis.py` - Complete workflow on a small dataset
- `02_conservation_analysis.py` - Detailed conservation scoring
- `03_motif_discovery.py` - Advanced motif discovery parameters

## Advanced Examples

- `04_phylogenetic_analysis.py` - Evolutionary analysis workflow
- `05_structural_mapping.py` - 3D structure integration
- `06_variant_analysis.py` - ClinVar variant enrichment

## Protein Family Examples

- `kinase_analysis.py` - Kinase domain motifs
- `gpcr_analysis.py` - GPCR conserved regions
- `transcription_factor.py` - DNA-binding motifs

## Running Examples

```bash
# Make sure EvoMotif is installed
pip install -e .

# Run an example
python examples/01_simple_analysis.py

# With custom parameters
python examples/02_conservation_analysis.py --protein TP53 --output results/
```

## Example Data

Small example datasets are provided in `example_data/`:

- `example_sequences.fasta` - 20 p53 sequences
- `example_alignment.fasta` - Pre-aligned sequences
- `example_structure.pdb` - Example PDB file
- `example_variants.csv` - Example variant data

## Notebooks

Jupyter notebooks with interactive examples:

- `tutorial_basic.ipynb` - Interactive basic tutorial
- `tutorial_advanced.ipynb` - Advanced features
- `visualization.ipynb` - Creating publication figures
