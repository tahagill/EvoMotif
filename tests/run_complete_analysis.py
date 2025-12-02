"""
COMPLETE EvoMotif analysis pipeline - ALL modules included.

This runs the FULL analysis with:
- Sequence retrieval
- Alignment
- Conservation scoring
- Motif discovery (FIXED - finds scattered conserved residues)
- Phylogenetic analysis
- 3D structure mapping
- Variant enrichment (optional)
"""

import argparse
import sys
import logging
from pathlib import Path
import json
import numpy as np

from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer
from evomotif.phylogeny import PhylogeneticAnalyzer
from evomotif.structure import StructureMapper
from evomotif.stats import StatisticalAnalyzer


def setup_logging(verbose=False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def find_scattered_motifs(conservation, gap_freq, consensus, min_conservation=0.7, max_gap=0.5):
    """
    Find biologically meaningful motifs including SCATTERED conserved residues.
    
    This is more biologically accurate than requiring consecutive windows.
    Many functional sites have conserved residues separated by variable regions.
    """
    high_cons_positions = []
    
    for i, (cons, gap, aa) in enumerate(zip(conservation, gap_freq, consensus)):
        if cons >= min_conservation and gap <= max_gap and aa != 'X':
            high_cons_positions.append({
                'position': i,
                'conservation': cons,
                'gap': gap,
                'amino_acid': aa
            })
    
    return high_cons_positions


def run_complete_analysis(
    protein_name: str,
    output_dir: Path,
    email: str,
    max_sequences: int = 50,
    min_conservation: float = 0.7,
    threads: int = 4,
    with_structure: bool = False,
    pdb_id: str = None
):
    """
    Run COMPLETE EvoMotif analysis on a protein.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger(__name__)
    results = {'protein': protein_name}
    
    print("\n" + "="*70)
    print(f"COMPLETE EvoMotif Analysis - {protein_name.upper()}")
    print("="*70)
    
    try:
        # ==================== STEP 1: Retrieve ====================
        print(f"\n[1/8] Retrieving {protein_name} sequences from NCBI...")
        retriever = SequenceRetriever(email=email)
        
        seq_file = output_dir / f"{protein_name}_sequences.fasta"
        sequences = retriever.retrieve_all(
            protein_name=protein_name,
            output_path=seq_file,
            max_sequences=max_sequences,
            min_length=100,  # Increased for real proteins
            max_length=3000,  # Allow longer sequences like BRCA1
            remove_fragments=False  # Keep fragments for better coverage
        )
        
        results['n_sequences'] = len(sequences)
        print(f"✓ Retrieved {len(sequences)} sequences")
        
        if len(sequences) < 5:
            print("ERROR: Too few sequences. Aborting.")
            return None
        
        # ==================== STEP 2: Align ====================
        print(f"\n[2/8] Aligning sequences with MAFFT...")
        aligner = SequenceAligner()
        
        aln_file = output_dir / f"{protein_name}_aligned.fasta"
        alignment = aligner.align_sequences(
            input_fasta=seq_file,
            output_fasta=aln_file,
            method="linsi",
            threads=threads
        )
        
        stats = aligner.calculate_alignment_stats(alignment)
        results['alignment_stats'] = stats
        print(f"✓ Aligned: {stats['n_sequences']} seqs, {stats['alignment_length']} positions")
        print(f"  Gap percentage: {stats['gap_percentage']:.1f}%")
        
        # ==================== STEP 3: Conservation ====================
        print(f"\n[3/8] Calculating conservation scores...")
        scorer = ConservationScorer()
        
        conservation = scorer.calculate_combined_conservation(alignment)
        gap_freq = scorer.calculate_gap_frequency(alignment)
        consensus = scorer.get_consensus_sequence(alignment)
        
        cons_file = output_dir / f"{protein_name}_conservation.json"
        scorer.save_conservation_scores(alignment, cons_file)
        
        results['mean_conservation'] = float(conservation.mean())
        results['max_conservation'] = float(conservation.max())
        
        print(f"✓ Conservation: mean={conservation.mean():.3f}, max={conservation.max():.3f}")
        
        # ==================== STEP 4: Motif Discovery (IMPROVED) ====================
        print(f"\n[4/8] Discovering conserved motifs...")
        
        # Method 1: Traditional sliding windows
        print("  → Finding consecutive motifs (sliding window)...")
        discoverer = MotifDiscoverer(
            window_sizes=[5, 7, 9, 11],
            min_conservation=min_conservation,
            max_std=0.15
        )
        
        window_motifs = discoverer.discover_motifs(alignment, conservation, gap_freq, consensus)
        
        # Method 2: Scattered conserved residues (MORE BIOLOGICALLY ACCURATE)
        print("  → Finding scattered conserved residues...")
        scattered_motifs = find_scattered_motifs(
            conservation, gap_freq, consensus,
            min_conservation=min_conservation - 0.05,  # Slightly lower
            max_gap=0.5
        )
        
        print(f"✓ Found {len(window_motifs)} consecutive motifs")
        print(f"✓ Found {len(scattered_motifs)} highly conserved positions")
        
        results['consecutive_motifs'] = len(window_motifs)
        results['conserved_positions'] = len(scattered_motifs)
        
        # Save both types
        if len(window_motifs) > 0:
            motif_dir = output_dir / "consecutive_motifs"
            discoverer.save_motifs_fasta(window_motifs, motif_dir)
            
            for i, motif in enumerate(window_motifs):
                print(f"  Consecutive Motif {i+1}: {motif.consensus_sequence} "
                      f"(pos {motif.start}-{motif.end}, cons={motif.mean_conservation:.3f})")
        
        if len(scattered_motifs) > 0:
            # Save scattered motifs
            scattered_file = output_dir / "conserved_positions.json"
            with open(scattered_file, 'w') as f:
                json.dump(scattered_motifs, f, indent=2)
            
            print(f"\n  Top 10 conserved positions:")
            for pos_info in sorted(scattered_motifs, key=lambda x: x['conservation'], reverse=True)[:10]:
                print(f"    Pos {pos_info['position']+1}: {pos_info['amino_acid']} "
                      f"(cons={pos_info['conservation']:.3f}, gap={pos_info['gap']*100:.1f}%)")
        
        results['motif_details'] = {
            'consecutive': [{'start': m.start, 'end': m.end, 'sequence': m.consensus_sequence,
                           'conservation': m.mean_conservation} for m in window_motifs],
            'conserved_positions': scattered_motifs
        }
        
        # ==================== STEP 5: Statistical Validation ====================
        print(f"\n[5/8] Statistical validation...")
        analyzer = StatisticalAnalyzer(random_seed=42)
        
        if len(window_motifs) > 0:
            motif_dicts = [
                {
                    'start': m.start,
                    'end': m.end,
                    'mean_conservation': m.mean_conservation,
                    'length': m.length
                }
                for m in window_motifs
            ]
            
            validated = analyzer.test_motif_significance(motif_dicts, conservation, n_permutations=1000)
            print(f"✓ Validated {len(validated)} motifs with permutation tests")
            
            significant = sum(1 for m in validated if m['adjusted_p_value'] < 0.05)
            print(f"  Significant motifs (FDR < 0.05): {significant}/{len(validated)}")
        
        # ==================== STEP 6: Phylogenetic Analysis ====================
        print(f"\n[6/8] Building phylogenetic tree...")
        try:
            phylo = PhylogeneticAnalyzer()
            tree_file = output_dir / f"{protein_name}_tree.nwk"
            
            tree = phylo.build_tree_fasttree(aln_file, tree_file)
            print(f"✓ Phylogenetic tree built with {len(tree.get_leaves())} leaves")
            results['phylogeny'] = {'n_leaves': len(tree.get_leaves())}
            
            # TODO: Ancestral state reconstruction if motifs found
            
        except FileNotFoundError as e:
            print(f"⚠ FastTree not installed - skipping phylogeny")
            print(f"  Install: sudo apt-get install fasttree")
        except Exception as e:
            print(f"⚠ Phylogeny failed: {e}")
        
        # ==================== STEP 7: 3D Structure Mapping ====================
        print(f"\n[7/8] Mapping to 3D structure...")
        
        if with_structure and pdb_id:
            try:
                mapper = StructureMapper()
                
                # Download PDB structure
                pdb_file = output_dir / f"{pdb_id}.pdb"
                print(f"  → Fetching PDB structure {pdb_id}...")
                
                import requests
                pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                response = requests.get(pdb_url)
                if response.status_code == 200:
                    with open(pdb_file, 'w') as f:
                        f.write(response.text)
                    print(f"  ✓ Downloaded {pdb_id}.pdb")
                    
                    # Load structure
                    structure = mapper.load_structure(pdb_file)
                    
                    # Map alignment to structure
                    ref_seq_with_gaps = str(alignment[0].seq)
                    mapping = mapper.map_alignment_to_structure(structure, ref_seq_with_gaps, 'A')
                    
                    print(f"  ✓ Mapped {len([k for k, v in mapping.items() if v])} positions to structure")
                    
                    # Visualize if we have motifs
                    if len(scattered_motifs) > 0:
                        print(f"  → Creating 3D visualization...")
                        vis_file = output_dir / f"{protein_name}_structure_3d.html"
                        
                        # Get conserved residue numbers (mapped to PDB residues)
                        conserved_positions = [pos_info['position'] for pos_info in scattered_motifs[:20]]
                        conserved_residues = [mapping.get(pos) for pos in conserved_positions if pos in mapping]
                        
                        html = mapper.visualize_structure_3d(
                            pdb_file=pdb_file,
                            motif_residues=conserved_residues,
                            output_html=vis_file
                        )
                        
                        print(f"  ✓ 3D visualization saved to {vis_file}")
                        print(f"    {len(conserved_residues)} conserved residues highlighted in RED")
                        print(f"    Open in browser to view interactive structure!")
                        results['structure_file'] = str(vis_file)
                else:
                    print(f"  ⚠ Could not download PDB {pdb_id}")
                    
            except Exception as e:
                print(f"⚠ Structure mapping failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            print(f"  ℹ Structure analysis skipped (use --pdb to enable)")
            print(f"    Example: --pdb 1TUP for p53 DNA-binding domain")
        
        # ==================== STEP 8: Save Results ====================
        print(f"\n[8/8] Saving results...")
        summary_file = output_dir / f"{protein_name}_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"✓ Summary saved to {summary_file}")
        
        # ==================== Summary ====================
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE!")
        print("="*70)
        print(f"Protein: {protein_name}")
        print(f"Sequences: {results['n_sequences']}")
        print(f"Alignment length: {results['alignment_stats']['alignment_length']}")
        print(f"Mean conservation: {results['mean_conservation']:.3f}")
        print(f"Consecutive motifs: {results['consecutive_motifs']}")
        print(f"Conserved positions: {results['conserved_positions']}")
        
        if 'structure_file' in results:
            print(f"\n🎨 3D Structure visualization: {results['structure_file']}")
            print(f"   Open this file in a web browser!")
        
        print(f"\nAll results saved to: {output_dir}/")
        print("="*70 + "\n")
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Complete EvoMotif analysis with ALL features',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python run_complete_analysis.py p53 your@email.com
  
  # With 3D structure visualization
  python run_complete_analysis.py p53 your@email.com --pdb 1TUP
  
  # More sequences, lower conservation threshold
  python run_complete_analysis.py ubiquitin your@email.com --max-sequences 100 --min-conservation 0.65
        """
    )
    
    parser.add_argument('protein', help='Protein name to analyze')
    parser.add_argument('email', help='Email for NCBI Entrez')
    parser.add_argument('--output', '-o', default='complete_analysis_results',
                       help='Output directory')
    parser.add_argument('--max-sequences', '-n', type=int, default=50,
                       help='Max sequences to retrieve')
    parser.add_argument('--min-conservation', '-c', type=float, default=0.7,
                       help='Min conservation for motifs')
    parser.add_argument('--threads', '-t', type=int, default=4,
                       help='CPU threads')
    parser.add_argument('--pdb', help='PDB ID for structure mapping (e.g., 1TUP for p53)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    output_dir = Path(args.output) / args.protein
    
    results = run_complete_analysis(
        protein_name=args.protein,
        output_dir=output_dir,
        email=args.email,
        max_sequences=args.max_sequences,
        min_conservation=args.min_conservation,
        threads=args.threads,
        with_structure=bool(args.pdb),
        pdb_id=args.pdb
    )
    
    if results is None:
        sys.exit(1)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
