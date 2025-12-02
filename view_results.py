#!/usr/bin/env python3
"""
Quick viewer for phylogenetic trees and conservation plots.
"""

import sys
import json
from pathlib import Path

def view_tree_ascii(tree_file):
    """Display tree in ASCII format."""
    try:
        from Bio import Phylo
        print(f"\n{'='*70}")
        print(f"Phylogenetic Tree: {tree_file.name}")
        print('='*70)
        tree = Phylo.read(tree_file, "newick")
        Phylo.draw_ascii(tree)
        print('='*70)
    except ImportError:
        print("⚠️  BioPython not installed. Install with: pip install biopython")
        print("\nShowing raw Newick format instead:")
        print(tree_file.read_text())

def view_tree_plot(tree_file, output_file=None):
    """Display tree as a plot."""
    try:
        from Bio import Phylo
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(12, 8))
        tree = Phylo.read(tree_file, "newick")
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title(f"Phylogenetic Tree: {tree_file.stem}", fontsize=16)
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✅ Tree plot saved to: {output_file}")
        else:
            plt.show()
            
    except ImportError as e:
        print(f"⚠️  Missing package: {e}")
        print("Install with: pip install biopython matplotlib")

def view_conservation_plot(json_file, output_file=None):
    """Plot conservation scores."""
    try:
        import matplotlib.pyplot as plt
        
        with open(json_file) as f:
            data = json.load(f)
        
        positions = [d['position'] for d in data]
        conservation = [d['conservation'] for d in data]
        gaps = [d['gap'] for d in data]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
        
        # Conservation plot
        ax1.plot(positions, conservation, 'b-', alpha=0.7, linewidth=1)
        ax1.fill_between(positions, conservation, alpha=0.3)
        ax1.axhline(y=0.7, color='r', linestyle='--', label='Threshold (0.7)', linewidth=2)
        ax1.set_ylabel('Conservation Score', fontsize=12)
        ax1.set_title(f'Conservation Profile: {json_file.parent.name}', fontsize=16, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1.05)
        
        # Gap frequency plot
        ax2.bar(positions, gaps, alpha=0.6, color='orange', width=1)
        ax2.set_xlabel('Position', fontsize=12)
        ax2.set_ylabel('Gap Frequency', fontsize=12)
        ax2.set_title('Gap Distribution', fontsize=14)
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.set_ylim(0, 1.05)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✅ Conservation plot saved to: {output_file}")
        else:
            plt.show()
            
    except ImportError as e:
        print(f"⚠️  Missing package: {e}")
        print("Install with: pip install matplotlib")

def show_summary(protein_dir):
    """Show analysis summary."""
    summary_file = protein_dir / f"{protein_dir.name}_summary.json"
    conserved_file = protein_dir / "conserved_positions.json"
    
    print(f"\n{'='*70}")
    print(f"📊 {protein_dir.name.upper()} ANALYSIS SUMMARY")
    print('='*70)
    
    if summary_file.exists():
        with open(summary_file) as f:
            summary = json.load(f)
        
        print(f"\n🧬 Sequences: {summary.get('n_sequences', 'N/A')}")
        
        if 'alignment_stats' in summary:
            print(f"📏 Alignment length: {summary['alignment_stats'].get('alignment_length', 'N/A')} bp")
        
        if 'conservation_stats' in summary:
            mean_cons = summary['conservation_stats'].get('mean', 'N/A')
            if isinstance(mean_cons, (int, float)):
                print(f"📊 Mean conservation: {mean_cons:.3f}")
            else:
                print(f"📊 Mean conservation: {mean_cons}")
        
        print(f"🔬 Conserved positions: {summary.get('n_conserved_positions', 'N/A')}")
        
        if 'motifs' in summary and 'n_consecutive' in summary['motifs']:
            print(f"🎯 Consecutive motifs: {summary['motifs']['n_consecutive']}")
    
    if conserved_file.exists():
        with open(conserved_file) as f:
            conserved = json.load(f)
        
        print(f"\n🌟 TOP 10 CONSERVED POSITIONS:")
        print(f"{'─'*70}")
        
        # Sort by conservation
        sorted_positions = sorted(conserved, key=lambda x: x['conservation'], reverse=True)[:10]
        
        for i, pos in enumerate(sorted_positions, 1):
            gap_pct = pos['gap'] * 100
            cons = pos['conservation']
            print(f"  {i:2d}. Pos {pos['position']:4d}: {pos['amino_acid']} "
                  f"(conservation={cons:.3f}, gap={gap_pct:4.1f}%)")
    
    print('='*70)

def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='View EvoMotif analysis results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # View AKT1 tree (ASCII)
  python view_results.py --tree complete_analysis_results/AKT1/
  
  # View P53 tree as plot
  python view_results.py --tree complete_analysis_results/p53/ --plot
  
  # View conservation profile
  python view_results.py --conservation complete_analysis_results/BRCA1/
  
  # Show summary
  python view_results.py --summary complete_analysis_results/AKT1/
  
  # Save plots to files
  python view_results.py --tree complete_analysis_results/p53/ --plot --output p53_tree.png
        """
    )
    
    parser.add_argument('protein_dir', nargs='?', 
                       help='Protein analysis directory (e.g., complete_analysis_results/AKT1/)')
    parser.add_argument('--tree', action='store_true', help='View phylogenetic tree')
    parser.add_argument('--conservation', action='store_true', help='View conservation plot')
    parser.add_argument('--summary', action='store_true', help='Show analysis summary')
    parser.add_argument('--plot', action='store_true', help='Show graphical plot (requires matplotlib)')
    parser.add_argument('--output', '-o', help='Save plot to file instead of displaying')
    parser.add_argument('--all', action='store_true', help='Show everything')
    
    args = parser.parse_args()
    
    # If no args, show help
    if not args.protein_dir:
        parser.print_help()
        print("\n📂 Available analyses:")
        results_dir = Path("complete_analysis_results")
        if results_dir.exists():
            for protein_dir in sorted(results_dir.iterdir()):
                if protein_dir.is_dir():
                    print(f"  - {protein_dir.name}")
        return
    
    protein_dir = Path(args.protein_dir)
    
    if not protein_dir.exists():
        print(f"❌ Error: Directory not found: {protein_dir}")
        return 1
    
    # Default to summary if nothing specified
    if not (args.tree or args.conservation or args.summary or args.all):
        args.summary = True
    
    # Show summary
    if args.summary or args.all:
        show_summary(protein_dir)
    
    # Show tree
    if args.tree or args.all:
        tree_file = list(protein_dir.glob("*.nwk"))
        if tree_file:
            tree_file = tree_file[0]
            if args.plot:
                output = args.output if args.output else None
                view_tree_plot(tree_file, output)
            else:
                view_tree_ascii(tree_file)
        else:
            print(f"⚠️  No tree file (.nwk) found in {protein_dir}")
    
    # Show conservation
    if args.conservation or args.all:
        conserved_file = protein_dir / "conserved_positions.json"
        if conserved_file.exists():
            output = args.output if args.output else None
            view_conservation_plot(conserved_file, output)
        else:
            print(f"⚠️  No conserved_positions.json found in {protein_dir}")

if __name__ == "__main__":
    sys.exit(main() or 0)
