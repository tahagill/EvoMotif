"""
Command-line interface for running EvoMotif integration tests.

This allows users to test any protein without writing code.
"""

import argparse
import sys
import logging
from pathlib import Path
import json

from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer
from evomotif.stats import StatisticalAnalyzer


def setup_logging(verbose=False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def run_analysis(
    protein_name: str,
    output_dir: Path,
    email: str,
    max_sequences: int = 50,
    min_conservation: float = 0.75,
    window_sizes: list = None,
    threads: int = 4
):
    """
    Run complete EvoMotif analysis on a protein.
    
    Args:
        protein_name: Name of protein to analyze
        output_dir: Directory to save results
        email: Email for NCBI Entrez
        max_sequences: Maximum sequences to retrieve
        min_conservation: Minimum conservation for motifs
        window_sizes: List of window sizes for motif discovery
        threads: Number of CPU threads
    
    Returns:
        Dictionary with analysis results
    """
    if window_sizes is None:
        window_sizes = [7, 9, 11]
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger(__name__)
    results = {}
    
    try:
        # Step 1: Retrieve sequences
        logger.info(f"Step 1/5: Retrieving {protein_name} sequences...")
        retriever = SequenceRetriever(email=email)
        
        seq_file = output_dir / f"{protein_name}_sequences.fasta"
        sequences = retriever.retrieve_all(
            protein_name=protein_name,
            output_path=seq_file,
            max_sequences=max_sequences,
            min_length=30,
            remove_fragments=True
        )
        
        results['n_sequences'] = len(sequences)
        logger.info(f"✓ Retrieved {len(sequences)} sequences")
        
        if len(sequences) < 5:
            logger.error("Too few sequences retrieved. Aborting.")
            return None
        
        # Step 2: Align sequences
        logger.info("Step 2/5: Aligning sequences...")
        aligner = SequenceAligner()
        
        aln_file = output_dir / f"{protein_name}_aligned.fasta"
        alignment = aligner.align_sequences(
            input_fasta=seq_file,
            output_fasta=aln_file,
            method="linsi",
            threads=threads
        )
        
        stats = aligner.calculate_alignment_stats(alignment)
        results['alignment_length'] = stats['alignment_length']
        results['gap_percentage'] = stats['gap_percentage']
        
        logger.info(f"✓ Alignment: {stats['n_sequences']} seqs, "
                   f"{stats['alignment_length']} positions, "
                   f"{stats['gap_percentage']:.1f}% gaps")
        
        # Step 3: Calculate conservation
        logger.info("Step 3/5: Calculating conservation...")
        scorer = ConservationScorer()
        
        conservation = scorer.calculate_combined_conservation(alignment)
        gap_freq = scorer.calculate_gap_frequency(alignment)
        consensus = scorer.get_consensus_sequence(alignment)
        
        cons_file = output_dir / f"{protein_name}_conservation.json"
        scorer.save_conservation_scores(alignment, cons_file)
        
        results['mean_conservation'] = float(conservation.mean())
        results['max_conservation'] = float(conservation.max())
        
        logger.info(f"✓ Conservation: mean={conservation.mean():.3f}, "
                   f"max={conservation.max():.3f}")
        
        # Step 4: Discover motifs
        logger.info("Step 4/5: Discovering motifs...")
        discoverer = MotifDiscoverer(
            window_sizes=window_sizes,
            min_conservation=min_conservation,
            max_std=0.1
        )
        
        motifs = discoverer.discover_motifs(
            alignment,
            conservation,
            gap_freq,
            consensus
        )
        
        results['n_motifs'] = len(motifs)
        
        if len(motifs) > 0:
            motif_dir = output_dir / "motifs"
            discoverer.save_motifs_fasta(motifs, motif_dir)
            
            logger.info(f"✓ Found {len(motifs)} motifs:")
            motif_details = []
            for i, motif in enumerate(motifs):
                logger.info(f"  Motif {i+1}: {motif.consensus_sequence} "
                          f"(pos {motif.start}-{motif.end}, "
                          f"cons={motif.mean_conservation:.3f})")
                
                motif_details.append({
                    'id': i + 1,
                    'sequence': motif.consensus_sequence,
                    'start': motif.start,
                    'end': motif.end,
                    'length': motif.length,
                    'conservation': motif.mean_conservation,
                    'n_sequences': len(motif.sequences)
                })
            
            results['motifs'] = motif_details
        else:
            logger.warning("No motifs found with current parameters")
            results['motifs'] = []
        
        # Step 5: Statistical validation
        logger.info("Step 5/5: Statistical validation...")
        analyzer = StatisticalAnalyzer(random_seed=42)
        
        if len(motifs) > 0:
            motif_dicts = [
                {
                    'start': m.start,
                    'end': m.end,
                    'mean_conservation': m.mean_conservation,
                    'length': m.length
                }
                for m in motifs
            ]
            
            validated = analyzer.test_motif_significance(
                motif_dicts,
                conservation,
                n_permutations=1000
            )
            
            logger.info("✓ Motif significance:")
            for i, motif in enumerate(validated):
                sig = "***" if motif['adjusted_p_value'] < 0.001 else \
                      "**" if motif['adjusted_p_value'] < 0.01 else \
                      "*" if motif['adjusted_p_value'] < 0.05 else "ns"
                
                logger.info(f"  Motif {i+1}: p={motif['p_value']:.4f}, "
                          f"adj_p={motif['adjusted_p_value']:.4f} {sig}")
                
                # Add significance to results
                results['motifs'][i]['p_value'] = motif['p_value']
                results['motifs'][i]['adjusted_p_value'] = motif['adjusted_p_value']
                results['motifs'][i]['z_score'] = motif['z_score']
        
        # Save summary
        summary_file = output_dir / f"{protein_name}_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"\n✅ Analysis complete!")
        logger.info(f"Results saved to: {output_dir}")
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=True)
        return None


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='EvoMotif: Test protein motif discovery pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze ubiquitin with default settings
  python run_integration_test.py ubiquitin myemail@example.com
  
  # Analyze p53 with custom parameters
  python run_integration_test.py p53 myemail@example.com \\
      --max-sequences 100 --min-conservation 0.8 \\
      --threads 8 --output results/p53
  
  # Quick test with fewer sequences
  python run_integration_test.py insulin myemail@example.com \\
      --max-sequences 20 --quick
        """
    )
    
    parser.add_argument(
        'protein',
        help='Name of protein to analyze (e.g., ubiquitin, p53, insulin)'
    )
    
    parser.add_argument(
        'email',
        help='Email address for NCBI Entrez (required)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='integration_test_results',
        help='Output directory (default: integration_test_results)'
    )
    
    parser.add_argument(
        '--max-sequences', '-n',
        type=int,
        default=50,
        help='Maximum sequences to retrieve (default: 50)'
    )
    
    parser.add_argument(
        '--min-conservation', '-c',
        type=float,
        default=0.75,
        help='Minimum conservation score for motifs (default: 0.75)'
    )
    
    parser.add_argument(
        '--window-sizes', '-w',
        type=int,
        nargs='+',
        default=[7, 9, 11],
        help='Window sizes for motif discovery (default: 7 9 11)'
    )
    
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=4,
        help='Number of CPU threads (default: 4)'
    )
    
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Quick mode: fewer sequences and permutations'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Adjust parameters for quick mode
    if args.quick:
        args.max_sequences = min(args.max_sequences, 20)
        logging.info("Quick mode: using max 20 sequences")
    
    # Create output directory with protein name
    output_dir = Path(args.output) / args.protein
    
    print("="*70)
    print(f"EvoMotif Integration Test - {args.protein.upper()}")
    print("="*70)
    print(f"Output directory: {output_dir}")
    print(f"Max sequences: {args.max_sequences}")
    print(f"Min conservation: {args.min_conservation}")
    print(f"Window sizes: {args.window_sizes}")
    print(f"Threads: {args.threads}")
    print("="*70 + "\n")
    
    # Run analysis
    results = run_analysis(
        protein_name=args.protein,
        output_dir=output_dir,
        email=args.email,
        max_sequences=args.max_sequences,
        min_conservation=args.min_conservation,
        window_sizes=args.window_sizes,
        threads=args.threads
    )
    
    if results is None:
        print("\n❌ Analysis failed. Check logs for details.")
        sys.exit(1)
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Protein: {args.protein}")
    print(f"Sequences retrieved: {results['n_sequences']}")
    print(f"Alignment length: {results['alignment_length']}")
    print(f"Mean conservation: {results['mean_conservation']:.3f}")
    print(f"Motifs found: {results['n_motifs']}")
    
    if results['n_motifs'] > 0:
        print(f"\nMotif Details:")
        for motif in results['motifs']:
            sig = "***" if motif.get('adjusted_p_value', 1) < 0.001 else \
                  "**" if motif.get('adjusted_p_value', 1) < 0.01 else \
                  "*" if motif.get('adjusted_p_value', 1) < 0.05 else "ns"
            
            print(f"  Motif {motif['id']}: {motif['sequence']}")
            print(f"    Position: {motif['start']}-{motif['end']}")
            print(f"    Conservation: {motif['conservation']:.3f}")
            if 'p_value' in motif:
                print(f"    P-value: {motif['p_value']:.4f} (adj: {motif['adjusted_p_value']:.4f}) {sig}")
    
    print(f"\nResults saved to: {output_dir}/")
    print("="*70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
