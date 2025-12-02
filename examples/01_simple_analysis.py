"""
Simple EvoMotif analysis example.

This script demonstrates a basic motif discovery workflow.
"""

import logging
from pathlib import Path
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def main():
    """Run simple motif discovery analysis."""
    
    # Create output directory
    output_dir = Path("output/simple_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("EvoMotif Simple Analysis Example")
    print("=" * 60)
    
    # Step 1: Retrieve sequences
    print("\n[1/5] Retrieving sequences...")
    retriever = SequenceRetriever(email="user@example.com")
    
    sequences = retriever.retrieve_all(
        protein_name="ubiquitin",
        output_path=output_dir / "sequences.fasta",
        max_sequences=50,
        min_length=50,
        remove_fragments=True
    )
    print(f"✓ Retrieved {len(sequences)} sequences")
    
    # Step 2: Align sequences
    print("\n[2/5] Aligning sequences...")
    aligner = SequenceAligner()
    
    alignment = aligner.align_sequences(
        input_fasta=output_dir / "sequences.fasta",
        output_fasta=output_dir / "alignment.fasta",
        method="linsi",
        threads=4
    )
    
    stats = aligner.calculate_alignment_stats(alignment)
    print(f"✓ Alignment: {stats['n_sequences']} seqs, {stats['alignment_length']} positions")
    
    # Step 3: Calculate conservation
    print("\n[3/5] Calculating conservation...")
    scorer = ConservationScorer()
    
    conservation = scorer.calculate_combined_conservation(alignment)
    gap_freq = scorer.calculate_gap_frequency(alignment)
    consensus = scorer.get_consensus_sequence(alignment)
    
    scorer.save_conservation_scores(
        alignment,
        output_dir / "conservation.json"
    )
    print(f"✓ Mean conservation: {conservation.mean():.3f}")
    
    # Step 4: Discover motifs
    print("\n[4/5] Discovering motifs...")
    discoverer = MotifDiscoverer(
        window_sizes=[7, 9, 11],
        min_conservation=0.85,
        max_std=0.1
    )
    
    motifs = discoverer.discover_motifs(
        alignment,
        conservation,
        gap_freq,
        consensus
    )
    
    print(f"✓ Found {len(motifs)} motifs:")
    for i, motif in enumerate(motifs):
        print(f"   Motif {i+1}: {motif.consensus_sequence} "
              f"(pos {motif.start}-{motif.end}, cons={motif.mean_conservation:.3f})")
    
    # Step 5: Save results
    print("\n[5/5] Saving results...")
    motif_dir = output_dir / "motifs"
    discoverer.save_motifs_fasta(motifs, motif_dir)
    
    print(f"\n✅ Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("\nGenerated files:")
    for file in sorted(output_dir.rglob("*")):
        if file.is_file():
            print(f"  - {file.relative_to(output_dir)}")


if __name__ == "__main__":
    main()
