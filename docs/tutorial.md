# Basic Tutorial: Discovering Motifs in p53

This tutorial demonstrates a complete EvoMotif workflow using the p53 tumor suppressor protein as an example.

## Overview

In this tutorial, you will:

1. Retrieve p53 sequences from multiple species
2. Perform multiple sequence alignment
3. Calculate conservation scores
4. Discover conserved motifs
5. Build HMM profiles
6. Analyze motif evolution
7. Map motifs to 3D structure
8. Analyze pathogenic variants

## Prerequisites

- EvoMotif installed with all dependencies
- Internet connection for sequence retrieval
- ~30 minutes to complete

## Step 1: Setup

```python
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

# Create output directory
output_dir = Path("tutorial_output")
output_dir.mkdir(exist_ok=True)
```

## Step 2: Retrieve Sequences

```python
from evomotif.retrieval import SequenceRetriever

# Initialize retriever (use your email)
retriever = SequenceRetriever(email="your.email@example.com")

# Define organisms of interest
organisms = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Danio rerio",
    "Xenopus laevis",
    "Drosophila melanogaster"
]

# Retrieve p53 sequences
sequences = retriever.retrieve_all(
    protein_name="p53",
    output_path=output_dir / "p53_sequences.fasta",
    organisms=organisms,
    max_sequences=100,
    min_length=200,
    remove_fragments=True
)

print(f"Retrieved {len(sequences)} sequences")
```

## Step 3: Multiple Sequence Alignment

```python
from evomotif.alignment import SequenceAligner

# Initialize aligner
aligner = SequenceAligner()

# Perform alignment using L-INS-i method
alignment = aligner.align_sequences(
    input_fasta=output_dir / "p53_sequences.fasta",
    output_fasta=output_dir / "p53_aligned.fasta",
    method="linsi",
    threads=4
)

# Calculate alignment statistics
stats = aligner.calculate_alignment_stats(alignment)
print(f"Alignment: {stats['n_sequences']} sequences, {stats['alignment_length']} positions")
print(f"Gap percentage: {stats['gap_percentage']:.2f}%")

# Optional: Trim poorly aligned regions
trimmed_alignment = aligner.trim_alignment(
    alignment,
    output_fasta=output_dir / "p53_trimmed.fasta",
    max_gap_fraction=0.5
)
```

## Step 4: Conservation Analysis

```python
from evomotif.conservation import ConservationScorer

# Initialize scorer
scorer = ConservationScorer()

# Calculate all conservation metrics
entropy = scorer.calculate_shannon_entropy(trimmed_alignment)
blosum = scorer.calculate_blosum_score(trimmed_alignment)
conservation = scorer.calculate_combined_conservation(trimmed_alignment)
gap_freq = scorer.calculate_gap_frequency(trimmed_alignment)
consensus = scorer.get_consensus_sequence(trimmed_alignment)

# Save conservation scores
scorer.save_conservation_scores(
    trimmed_alignment,
    output_dir / "conservation_scores.json"
)

# Visualize conservation
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 4))
plt.plot(conservation, label='Combined Conservation')
plt.axhline(y=0.85, color='r', linestyle='--', label='Motif Threshold')
plt.xlabel('Alignment Position')
plt.ylabel('Conservation Score')
plt.title('p53 Conservation Profile')
plt.legend()
plt.savefig(output_dir / "conservation_plot.png", dpi=300)
print("Conservation plot saved")
```

## Step 5: Motif Discovery

```python
from evomotif.motif_discovery import MotifDiscoverer, HMMProfileBuilder

# Initialize discoverer
discoverer = MotifDiscoverer(
    window_sizes=[7, 9, 11, 13],
    min_conservation=0.85,
    max_std=0.1
)

# Discover motifs
motifs = discoverer.discover_motifs(
    trimmed_alignment,
    conservation,
    gap_freq,
    consensus
)

print(f"\nDiscovered {len(motifs)} motifs:")
for i, motif in enumerate(motifs):
    print(f"  Motif {i+1}: position {motif.start}-{motif.end}, "
          f"conservation={motif.mean_conservation:.3f}, "
          f"sequence={motif.consensus_sequence}")

# Save motif sequences
motif_dir = output_dir / "motifs"
discoverer.save_motifs_fasta(motifs, motif_dir)

# Build HMM profiles
hmm_builder = HMMProfileBuilder()
hmm_dir = output_dir / "hmm_profiles"
motifs = hmm_builder.build_all_profiles(motifs, motif_dir, hmm_dir)

# Search for motifs in original sequences
for i, motif in enumerate(motifs):
    if motif.hmm_path:
        hits = hmm_builder.search_sequences(
            motif.hmm_path,
            output_dir / "p53_sequences.fasta",
            output_dir / f"motif_{i+1}_hits.txt"
        )
        print(f"Motif {i+1}: {len(hits)} hits")
```

## Step 6: Statistical Validation

```python
from evomotif.stats import StatisticalAnalyzer

# Initialize analyzer
analyzer = StatisticalAnalyzer()

# Test motif significance vs random
motif_dicts = [
    {
        'start': m.start,
        'end': m.end,
        'mean_conservation': m.mean_conservation,
        'length': m.length
    }
    for m in motifs
]

significant_motifs = analyzer.test_motif_significance(
    motif_dicts,
    conservation,
    n_permutations=1000
)

print("\nMotif significance:")
for i, motif in enumerate(significant_motifs):
    sig = "***" if motif['adjusted_p_value'] < 0.001 else \
          "**" if motif['adjusted_p_value'] < 0.01 else \
          "*" if motif['adjusted_p_value'] < 0.05 else "ns"
    print(f"  Motif {i+1}: p={motif['adjusted_p_value']:.4f} {sig}")
```

## Step 7: Phylogenetic Analysis

```python
from evomotif.phylogeny import PhylogeneticAnalyzer

# Initialize analyzer
phylo = PhylogeneticAnalyzer()

# Build phylogenetic tree
tree = phylo.build_tree_fasttree(
    output_dir / "p53_trimmed.fasta",
    output_dir / "p53_tree.nwk"
)

# Create motif presence/absence matrix
motif_presence = {}
for record in trimmed_alignment:
    motif_presence[record.id] = {
        i: 1  # Simplification: all sequences have all motifs
        for i in range(len(motifs))
    }

# Annotate tree with motif presence
tree = phylo.annotate_tree_with_motifs(tree, motif_presence)

# Reconstruct ancestral states for each motif
for i in range(len(motifs)):
    tree = phylo.reconstruct_ancestral_states(tree, i)
    origin = phylo.identify_motif_origin(tree, i)
    print(f"Motif {i+1} originated at: {origin}")

# Visualize tree
phylo.visualize_tree(
    tree,
    output_dir / "phylo_tree.png",
    motif_id=0
)
```

## Step 8: Structural Mapping

```python
from evomotif.structure import StructureMapper

# Initialize mapper
mapper = StructureMapper()

# Download p53 structure (PDB: 1TUP)
# Assume you have the PDB file
pdb_file = Path("data/1tup.pdb")

if pdb_file.exists():
    # Load structure
    structure = mapper.load_structure(pdb_file)
    
    # Map alignment to structure
    reference_seq = str(trimmed_alignment[0].seq)
    aln_to_struct = mapper.map_alignment_to_structure(
        structure,
        reference_seq,
        chain_id='A'
    )
    
    # Map conservation to structure
    cons_map = mapper.map_conservation_to_structure(
        conservation,
        aln_to_struct
    )
    
    # Calculate residue properties
    residue_props = mapper.calculate_residue_properties(
        structure,
        chain_id='A'
    )
    
    # Identify motif clusters in 3D
    motif_residues = []
    for motif in motifs:
        for pos in range(motif.start, motif.end):
            if pos in aln_to_struct:
                motif_residues.append(aln_to_struct[pos])
    
    clusters = mapper.identify_motif_clusters(
        motif_residues,
        structure,
        chain_id='A',
        distance_threshold=8.0
    )
    
    print(f"\nIdentified {len(clusters)} spatial clusters")
    
    # Test clustering significance
    clustering_stats = mapper.calculate_clustering_statistics(
        motif_residues,
        structure,
        chain_id='A',
        n_permutations=1000
    )
    
    print(f"Clustering p-value: {clustering_stats['p_value']:.4f}")
    
    # Create 3D visualization
    mapper.visualize_structure_3d(
        pdb_file,
        conservation_map=cons_map,
        motif_residues=motif_residues,
        output_html=output_dir / "structure_3d.html"
    )
    
    print("3D visualization saved to structure_3d.html")
```

## Step 9: Variant Analysis

```python
from evomotif.variants import VariantAnalyzer
import pandas as pd

# Initialize analyzer
var_analyzer = VariantAnalyzer()

# Load ClinVar variants (example data)
# You would download this from ClinVar
variants_example = pd.DataFrame({
    'position': [175, 245, 248, 273, 282],
    'clinical_significance': ['pathogenic', 'pathogenic', 'pathogenic', 'pathogenic', 'benign'],
    'variant': ['R175H', 'G245S', 'R248Q', 'R273H', 'R282W']
})

# Map variants to motifs
motif_positions = [(m.start, m.end) for m in motifs]
variants_mapped = var_analyzer.map_variants_to_motifs(
    variants_example,
    motif_positions
)

# Test enrichment
enrichment = var_analyzer.test_variant_enrichment(
    variants_mapped,
    total_positions=len(consensus),
    motif_length=sum(m.length for m in motifs)
)

print("\nVariant enrichment analysis:")
print(f"  Odds ratio: {enrichment['odds_ratio']:.2f}")
print(f"  P-value: {enrichment['p_value']:.4f}")
print(f"  Fold enrichment: {enrichment['enrichment_fold']:.2f}x")

# Permutation test
perm_result = var_analyzer.permutation_test_enrichment(
    variants_mapped,
    motif_positions,
    total_positions=len(consensus),
    n_permutations=10000
)

print(f"  Permutation p-value: {perm_result['p_value']:.4f}")

# Per-motif analysis
motif_variants = var_analyzer.analyze_motif_specific_variants(
    variants_mapped,
    motif_positions
)

for mv in motif_variants:
    print(f"  Motif {mv['motif_id']+1}: {mv['pathogenic']} pathogenic variants")
```

## Step 10: Generate Report

```python
# Save comprehensive analysis
from evomotif.stats import StatisticalAnalyzer

analyzer.save_statistical_report(
    output_dir / "final_report.json",
    motifs=significant_motifs,
    enrichment=enrichment,
    permutation=perm_result,
    clustering=clustering_stats if pdb_file.exists() else {}
)

print(f"\n✅ Analysis complete! Results saved to {output_dir}/")
print("\nGenerated files:")
for file in sorted(output_dir.rglob("*")):
    if file.is_file():
        print(f"  - {file.relative_to(output_dir)}")
```

## Expected Results

For p53, you should find:

1. **3-5 highly conserved motifs**, including:
   - DNA-binding domain residues (positions ~100-300)
   - Zinc-binding motif
   - Oligomerization domain

2. **Significant conservation** (mean > 0.85) in functional regions

3. **Spatial clustering** of motifs in 3D structure (p < 0.05)

4. **Variant enrichment** in motifs (OR > 2, p < 0.01)

5. **Evolutionary origin** in early eukaryotes

## Next Steps

- Try different proteins (kinases, GPCRs, transcription factors)
- Adjust motif discovery parameters
- Compare with Pfam domains
- Analyze motif co-evolution
- Integrate with AlphaFold structures

## Troubleshooting

**Issue**: Alignment takes too long
- Use `method="auto"` instead of `"linsi"`
- Reduce number of sequences with CD-HIT

**Issue**: No motifs found
- Lower `min_conservation` threshold (try 0.75-0.80)
- Increase `max_std` (try 0.15-0.20)
- Check alignment quality

**Issue**: Memory error
- Reduce number of sequences
- Use trimmed alignment
- Close other applications

## Resources

- [API Documentation](../api.md)
- [Advanced Examples](../examples/index.md)
- [GitHub Issues](https://github.com/yourusername/EvoMotif/issues)
