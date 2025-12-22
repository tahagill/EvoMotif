"""
High-level pipeline interface for EvoMotif.

Provides simple, scientist-friendly API for complete protein analysis
without requiring detailed knowledge of individual modules.
"""

import logging
import shutil
import sys
from pathlib import Path
from typing import Optional, Dict, List, Any
import numpy as np
import json

from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer
from evomotif.stats import StatisticalAnalyzer
from evomotif.phylogeny import PhylogeneticAnalyzer
from evomotif.structure import StructureMapper

logger = logging.getLogger(__name__)


class AnalysisResults:
    """Container for analysis results with convenient access methods."""
    
    def __init__(self, data: Dict[str, Any]):
        """Initialize results container."""
        self.data = data
        self.protein = data.get('protein', 'Unknown')
        self.n_sequences = data.get('n_sequences', 0)
        self.motifs = data.get('motifs', [])
        self.conserved_positions = data.get('conserved_positions', [])
        self.conservation_scores = data.get('conservation', None)
        self.output_dir = Path(data.get('output_dir', '.'))
        self.files = data.get('files', {})
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "="*60,
            f"EvoMotif Analysis Results: {self.protein}",
            "="*60,
            f"Sequences analyzed: {self.n_sequences}",
            f"Consecutive motifs found: {len(self.motifs)}",
            f"Conserved positions: {len(self.conserved_positions)}",
        ]
        
        if self.conservation_scores is not None:
            lines.append(f"Mean conservation: {np.mean(self.conservation_scores):.3f}")
            lines.append(f"Max conservation: {np.max(self.conservation_scores):.3f}")
        
        if self.motifs:
            lines.append("\nTop 5 motifs:")
            for i, motif in enumerate(self.motifs[:5], 1):
                lines.append(f"  {i}. {motif.get('sequence', 'N/A')} "
                           f"(pos {motif.get('start', '?')}-{motif.get('end', '?')}, "
                           f"cons={motif.get('conservation', 0):.3f})")
        
        lines.append(f"\n📁 Results saved to: {self.output_dir}")
        lines.append("="*60)
        
        return "\n".join(lines)
    
    def get_file(self, file_type: str) -> Optional[Path]:
        """Get path to a specific output file."""
        return self.files.get(file_type)
    
    def export_json(self, filepath: Path) -> None:
        """Export results to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.data, f, indent=2, default=str)
    
    def __repr__(self):
        return f"AnalysisResults(protein={self.protein}, motifs={len(self.motifs)}, sequences={self.n_sequences})"


class EvoMotifPipeline:
    """
    High-level interface for complete protein motif analysis.
    
    This class wraps all EvoMotif modules into a simple, easy-to-use API
    that requires minimal configuration while providing sensible defaults.
    """
    
    def __init__(self):
        """Initialize pipeline."""
        self.logger = logging.getLogger(__name__)
    
    @staticmethod
    def check_dependencies(required: List[str] = None) -> Dict[str, bool]:
        """
        Check if required external tools are installed.
        
        Args:
            required: List of tools to check. Default: ['mafft', 'fasttree']
        
        Returns:
            Dictionary mapping tool names to availability (True/False)
        """
        if required is None:
            required = ['mafft', 'fasttree']
        
        status = {}
        for tool in required:
            status[tool] = shutil.which(tool) is not None
        
        return status
    
    def _check_and_report_dependencies(self, include_structure: bool = False) -> bool:
        """Check dependencies and print helpful error messages if missing."""
        tools = ['mafft', 'fasttree']
        if include_structure:
            tools.append('dssp')
        
        status = self.check_dependencies(tools)
        missing = [tool for tool, available in status.items() if not available]
        
        if missing:
            print("\n⚠️  Missing required tools:")
            for tool in missing:
                print(f"  ❌ {tool}")
            
            print("\n📦 Installation instructions:")
            print("  Ubuntu/Debian:")
            print(f"    sudo apt-get install {' '.join(missing)}")
            print("  macOS:")
            print(f"    brew install {' '.join(missing)}")
            print("  Conda:")
            print(f"    conda install -c bioconda {' '.join(missing)}")
            
            return False
        
        return True
    
    def analyze(
        self,
        protein_name: str,
        email: str,
        output_dir: Optional[str] = None,
        pdb_id: Optional[str] = None,
        max_sequences: int = 50,
        min_conservation: float = 0.70,
        threads: int = 4,
        verbose: bool = False,
        check_deps: bool = True
    ) -> AnalysisResults:
        """
        Run complete EvoMotif analysis with one function call.
        
        This is the main entry point for most users. It runs the complete
        pipeline from sequence retrieval to final analysis and returns
        an easy-to-use results object.
        
        Args:
            protein_name: Name of protein to analyze (e.g., "p53", "ubiquitin")
            email: Email for NCBI Entrez (required by NCBI)
            output_dir: Output directory path (default: "./evomotif_results/{protein}")
            pdb_id: Optional PDB ID for structure mapping (e.g., "1TUP")
            max_sequences: Maximum sequences to retrieve (default: 50)
            min_conservation: Conservation threshold 0-1 (default: 0.70)
            threads: CPU threads for alignment (default: 4)
            verbose: Enable debug logging (default: False)
            check_deps: Check for required tools before running (default: True)
        
        Returns:
            AnalysisResults object with all results and file paths
        
        Raises:
            RuntimeError: If required dependencies are missing
            ValueError: If parameters are invalid
        
        Examples:
            >>> import evomotif
            >>> results = evomotif.analyze_protein("p53", "user@email.com")
            >>> print(results.summary())
            >>> 
            >>> # With structure
            >>> results = evomotif.analyze_protein("p53", "user@email.com", pdb_id="1TUP")
            >>> 
            >>> # Custom parameters
            >>> results = evomotif.analyze_protein(
            ...     "BRCA1", "user@email.com",
            ...     max_sequences=100,
            ...     min_conservation=0.65,
            ...     threads=8
            ... )
        """
        # Setup logging
        level = logging.DEBUG if verbose else logging.INFO
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        
        # Validate parameters
        if not protein_name:
            raise ValueError("protein_name is required")
        if not email:
            raise ValueError("email is required for NCBI Entrez")
        if not 0 <= min_conservation <= 1:
            raise ValueError("min_conservation must be between 0 and 1")
        if max_sequences < 5:
            raise ValueError("max_sequences must be at least 5")
        
        # Check dependencies
        if check_deps:
            if not self._check_and_report_dependencies(include_structure=bool(pdb_id)):
                raise RuntimeError(
                    "Required dependencies are missing. "
                    "Please install them before running analysis."
                )
        
        # Setup output directory
        if output_dir is None:
            output_dir = f"./evomotif_results/{protein_name}"
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results_data = {
            'protein': protein_name,
            'output_dir': str(output_dir),
            'files': {}
        }
        
        print("\n" + "="*70)
        print(f"🧬 EvoMotif Analysis: {protein_name.upper()}")
        print("="*70)
        
        try:
            # Step 1: Retrieve sequences
            print(f"\n[1/6] 🔍 Retrieving sequences from NCBI...")
            retriever = SequenceRetriever(email=email)
            
            seq_file = output_dir / f"{protein_name}_sequences.fasta"
            sequences = retriever.retrieve_all(
                protein_name=protein_name,
                output_path=seq_file,
                max_sequences=max_sequences,
                min_length=100,
                max_length=5000,
                remove_fragments=False
            )
            
            results_data['n_sequences'] = len(sequences)
            results_data['files']['sequences'] = seq_file
            print(f"      ✓ Retrieved {len(sequences)} sequences")
            
            if len(sequences) < 5:
                raise RuntimeError(
                    f"Too few sequences ({len(sequences)}). "
                    "Try increasing max_sequences or using a different protein."
                )
            
            # Step 2: Align sequences
            print(f"\n[2/6] 🧩 Aligning sequences with MAFFT...")
            aligner = SequenceAligner()
            
            aln_file = output_dir / f"{protein_name}_aligned.fasta"
            alignment = aligner.align_sequences(
                input_fasta=seq_file,
                output_fasta=aln_file,
                method="linsi",
                threads=threads
            )
            
            stats = aligner.calculate_alignment_stats(alignment)
            results_data['alignment_stats'] = stats
            results_data['files']['alignment'] = aln_file
            print(f"      ✓ Aligned: {stats['n_sequences']} sequences, "
                  f"{stats['alignment_length']} positions "
                  f"({stats['gap_percentage']:.1f}% gaps)")
            
            # Step 3: Calculate conservation
            print(f"\n[3/6] 📊 Calculating conservation scores...")
            scorer = ConservationScorer()
            
            conservation = scorer.calculate_combined_conservation(alignment)
            gap_freq = scorer.calculate_gap_frequency(alignment)
            consensus = scorer.get_consensus_sequence(alignment)
            
            cons_file = output_dir / f"{protein_name}_conservation.json"
            scorer.save_conservation_scores(alignment, cons_file)
            
            results_data['conservation'] = conservation.tolist()
            results_data['mean_conservation'] = float(conservation.mean())
            results_data['max_conservation'] = float(conservation.max())
            results_data['files']['conservation'] = cons_file
            print(f"      ✓ Mean conservation: {conservation.mean():.3f}, "
                  f"Max: {conservation.max():.3f}")
            
            # Step 4: Discover motifs
            print(f"\n[4/6] 🎯 Discovering conserved motifs...")
            
            # Window-based motifs
            discoverer = MotifDiscoverer(
                window_sizes=[5, 7, 9, 11],
                min_conservation=min_conservation,
                max_std=0.15
            )
            
            motifs = discoverer.discover_motifs(alignment, conservation, gap_freq, consensus)
            
            # Scattered conserved positions
            scattered = []
            for i, (cons, gap, aa) in enumerate(zip(conservation, gap_freq, consensus)):
                if cons >= min_conservation - 0.05 and gap <= 0.5 and aa != 'X':
                    scattered.append({
                        'position': int(i),
                        'conservation': float(cons),
                        'gap_frequency': float(gap),
                        'amino_acid': aa
                    })
            
            results_data['motifs'] = [
                {
                    'start': m.start,
                    'end': m.end,
                    'sequence': m.consensus_sequence,
                    'conservation': m.mean_conservation
                }
                for m in motifs
            ]
            results_data['conserved_positions'] = scattered
            
            print(f"      ✓ Found {len(motifs)} consecutive motifs")
            print(f"      ✓ Found {len(scattered)} highly conserved positions")
            
            if motifs:
                motif_dir = output_dir / "motifs"
                discoverer.save_motifs_fasta(motifs, motif_dir)
                results_data['files']['motifs_dir'] = motif_dir
            
            if scattered:
                scattered_file = output_dir / "conserved_positions.json"
                with open(scattered_file, 'w') as f:
                    json.dump(scattered, f, indent=2)
                results_data['files']['conserved_positions'] = scattered_file
            
            # Step 5: Statistical validation
            print(f"\n[5/6] 📈 Statistical validation...")
            analyzer = StatisticalAnalyzer(random_seed=42)
            
            if motifs:
                motif_dicts = [
                    {
                        'start': m.start,
                        'end': m.end,
                        'mean_conservation': m.mean_conservation
                    }
                    for m in motifs
                ]
                
                validated = analyzer.test_motif_significance(
                    motif_dicts,
                    conservation,
                    n_permutations=1000
                )
                
                results_data['validated_motifs'] = validated
                n_significant = sum(1 for m in validated if m['adjusted_p_value'] < 0.05)
                print(f"      ✓ {n_significant}/{len(validated)} motifs statistically significant (FDR < 0.05)")
            
            # Step 6: Phylogenetic analysis
            print(f"\n[6/6] 🌳 Building phylogenetic tree...")
            phylo = PhylogeneticAnalyzer()
            
            tree_file = output_dir / f"{protein_name}_tree.nwk"
            try:
                tree = phylo.build_tree_fasttree(
                    alignment_file=aln_file,
                    output_tree=tree_file,
                    model="JTT"
                )
                results_data['files']['tree'] = tree_file
                print(f"      ✓ Phylogenetic tree saved")
            except Exception as e:
                print(f"      ⚠️  Tree building failed: {e}")
            
            # Optional: Structure mapping
            if pdb_id:
                print(f"\n[BONUS] 🧊 Mapping to 3D structure (PDB: {pdb_id})...")
                try:
                    mapper = StructureMapper()
                    from Bio.PDB import PDBList
                    
                    pdb_file = output_dir / f"{pdb_id}.pdb"
                    if not pdb_file.exists():
                        pdbl = PDBList()
                        pdbl.retrieve_pdb_file(pdb_id, pdir=str(output_dir), file_format='pdb')
                        # Rename downloaded file
                        downloaded = output_dir / f"pdb{pdb_id.lower()}.ent"
                        if downloaded.exists():
                            downloaded.rename(pdb_file)
                    
                    structure = mapper.load_structure(pdb_file)
                    results_data['files']['pdb'] = pdb_file
                    print(f"          ✓ Structure mapped successfully")
                    
                except Exception as e:
                    print(f"          ⚠️  Structure mapping failed: {e}")
            
            # Save summary
            summary_file = output_dir / f"{protein_name}_summary.json"
            with open(summary_file, 'w') as f:
                json.dump(results_data, f, indent=2, default=str)
            results_data['files']['summary'] = summary_file
            
            print("\n" + "="*70)
            print("✅ Analysis complete!")
            print(f"📁 Results saved to: {output_dir}")
            print("="*70 + "\n")
            
            return AnalysisResults(results_data)
            
        except Exception as e:
            self.logger.error(f"Analysis failed: {e}", exc_info=verbose)
            raise RuntimeError(f"Analysis failed: {e}") from e


# Convenience function for direct import
def analyze_protein(*args, **kwargs) -> AnalysisResults:
    """
    Run complete protein motif analysis (convenience function).
    
    This is a shortcut for EvoMotifPipeline().analyze().
    See EvoMotifPipeline.analyze() for full documentation.
    
    Examples:
        >>> import evomotif
        >>> results = evomotif.analyze_protein("p53", "user@email.com")
        >>> print(results.summary())
    """
    pipeline = EvoMotifPipeline()
    return pipeline.analyze(*args, **kwargs)
