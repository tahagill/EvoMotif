"""
Structural mapping module.

Maps motifs and conservation scores onto 3D protein structures
using PDB files or AlphaFold models.
"""

import numpy as np
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select, DSSP
from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import aa1, aa3, d3_to_index, dindex_to_1
import py3Dmol
import json

# Create 3-letter to 1-letter amino acid mapping
THREE_TO_ONE = {three: dindex_to_1.get(idx, 'X') for three, idx in d3_to_index.items()}

logger = logging.getLogger(__name__)


class StructureMapper:
    """Map conservation scores and motifs onto 3D structures."""
    
    def __init__(self):
        """Initialize structure mapper."""
        self.logger = logging.getLogger(__name__)
        self.parser = PDBParser(QUIET=True)
        self.pdb_file_path = None  # Store for DSSP
    
    def load_structure(
        self,
        pdb_file: Path,
        structure_id: Optional[str] = None
    ) -> Structure:
        """
        Load PDB structure file.
        
        Args:
            pdb_file: Path to PDB file
            structure_id: Optional structure identifier
            
        Returns:
            BioPython Structure object
        """
        if structure_id is None:
            structure_id = pdb_file.stem
        
        self.pdb_file_path = pdb_file  # Store for DSSP
        structure = self.parser.get_structure(structure_id, pdb_file)
        self.logger.info(f"Loaded structure: {structure_id}")
        
        return structure
    
    def get_alphafold_confidence(
        self,
        structure: Structure,
        chain_id: str = 'A'
    ) -> Dict[int, float]:
        """
        Extract pLDDT confidence scores from AlphaFold model.
        
        AlphaFold models store per-residue confidence (pLDDT) in the 
        B-factor field of the PDB file. This method extracts those values.
        
        Args:
            structure: BioPython Structure object (AlphaFold model)
            chain_id: PDB chain identifier
            
        Returns:
            Dictionary mapping {residue_number: pLDDT_score}
            
        Note:
            pLDDT scores range from 0-100:
            - >90: Very high confidence
            - 70-90: Confident
            - 50-70: Low confidence
            - <50: Very low confidence (disordered)
        """
        confidence = {}
        
        try:
            chain = structure[0][chain_id]
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    # Get B-factor from CA atom (pLDDT in AlphaFold models)
                    if 'CA' in residue:
                        plddt = residue['CA'].get_bfactor()
                        confidence[residue.id[1]] = float(plddt)
            
            self.logger.info(f"Extracted pLDDT scores for {len(confidence)} residues")
        except Exception as e:
            self.logger.error(f"Failed to extract AlphaFold confidence: {e}")
        
        return confidence
    
    def map_alignment_to_structure(
        self,
        structure: Structure,
        alignment_seq: str,
        chain_id: str = 'A'
    ) -> Dict[int, int]:
        """
        Map alignment positions to structure residue numbers.
        
        Uses positional mapping. For robust mapping with indels,
        consider using pairwise alignment (e.g., Bio.pairwise2).
        
        Args:
            structure: BioPython Structure object
            alignment_seq: Sequence from alignment (with gaps)
            chain_id: PDB chain identifier
            
        Returns:
            Dictionary mapping {alignment_pos: residue_number}
        """
        # Get structure sequence
        chain = structure[0][chain_id]
        structure_residues = [res for res in chain if res.id[0] == ' ']
        
        structure_seq = ''.join([
            THREE_TO_ONE.get(res.resname, 'X') for res in structure_residues
        ])
        
        # Remove gaps from alignment sequence
        ungapped_seq = alignment_seq.replace('-', '')
        
        # Simple alignment (should use proper alignment for real case)
        if ungapped_seq not in structure_seq:
            self.logger.warning(
                "Alignment sequence not found in structure. "
                "Using positional mapping."
            )
        
        # Create mapping
        mapping = {}
        aln_pos = 0
        struct_pos = 0
        
        for i, aa in enumerate(alignment_seq):
            if aa != '-':
                if struct_pos < len(structure_residues):
                    res_id = structure_residues[struct_pos].id[1]
                    mapping[i] = res_id
                    struct_pos += 1
            aln_pos += 1
        
        self.logger.info(f"Mapped {len(mapping)} alignment positions to structure")
        return mapping
    
    def calculate_residue_properties(
        self,
        structure: Structure,
        chain_id: str = 'A',
        dssp_path: str = 'dssp'
    ) -> Dict[int, Dict]:
        """
        Calculate structural properties for each residue.
        
        Args:
            structure: BioPython Structure object
            chain_id: PDB chain identifier
            dssp_path: Path to DSSP executable
            
        Returns:
            Dictionary mapping {residue_number: properties}
        """
        properties = {}
        chain = structure[0][chain_id]
        
        # Try to run DSSP for secondary structure and accessibility
        dssp_data = None
        try:
            if self.pdb_file_path:
                dssp_data = DSSP(structure[0], str(self.pdb_file_path), dssp=dssp_path)
                self.logger.info("DSSP analysis completed")
            else:
                self.logger.warning("PDB file path not stored; cannot run DSSP")
        except Exception as e:
            self.logger.warning(f"DSSP analysis failed: {e}")
        
        for residue in chain:
            if residue.id[0] != ' ':  # Skip hetero atoms
                continue
            
            res_num = residue.id[1]
            res_props = {
                'residue': THREE_TO_ONE.get(residue.resname, 'X'),
                'secondary_structure': None,
                'accessibility': None,
                'coordinates': None
            }
            
            # Get DSSP data if available
            if dssp_data:
                try:
                    dssp_key = (chain_id, residue.id)
                    if dssp_key in dssp_data:
                        dssp_res = dssp_data[dssp_key]
                        res_props['secondary_structure'] = dssp_res[2]  # H, E, C
                        res_props['accessibility'] = dssp_res[3]  # RSA
                except:
                    pass
            
            # Get CA coordinates
            if 'CA' in residue:
                ca_atom = residue['CA']
                res_props['coordinates'] = ca_atom.get_coord().tolist()
            
            properties[res_num] = res_props
        
        return properties
    
    def map_conservation_to_structure(
        self,
        conservation: np.ndarray,
        alignment_to_structure: Dict[int, int]
    ) -> Dict[int, float]:
        """
        Map conservation scores to structure residues.
        
        Args:
            conservation: Conservation score array
            alignment_to_structure: Mapping from alignment to structure
            
        Returns:
            Dictionary mapping {residue_number: conservation_score}
        """
        conservation_map = {}
        
        for aln_pos, struct_res in alignment_to_structure.items():
            if aln_pos < len(conservation):
                conservation_map[struct_res] = float(conservation[aln_pos])
        
        self.logger.info(
            f"Mapped conservation scores to {len(conservation_map)} residues"
        )
        return conservation_map
    
    def identify_motif_clusters(
        self,
        motif_residues: List[int],
        structure: Structure,
        chain_id: str = 'A',
        distance_threshold: float = 8.0
    ) -> List[List[int]]:
        """
        Identify spatial clusters of motif residues.
        
        Args:
            motif_residues: List of motif residue numbers
            structure: BioPython Structure object
            chain_id: PDB chain identifier
            distance_threshold: Maximum CA-CA distance for clustering (Angstroms)
            
        Returns:
            List of residue clusters
        """
        chain = structure[0][chain_id]
        
        # Get CA coordinates for motif residues
        coords = {}
        for res_num in motif_residues:
            try:
                residue = chain[res_num]
                if 'CA' in residue:
                    coords[res_num] = residue['CA'].get_coord()
            except KeyError:
                continue
        
        if not coords:
            return []
        
        # Simple clustering using distance threshold
        clusters = []
        unclustered = set(coords.keys())
        
        while unclustered:
            # Start new cluster with arbitrary residue
            seed = unclustered.pop()
            cluster = [seed]
            
            # Add nearby residues
            changed = True
            while changed:
                changed = False
                for res_num in list(unclustered):
                    # Check distance to any cluster member
                    for cluster_res in cluster:
                        dist = np.linalg.norm(
                            coords[res_num] - coords[cluster_res]
                        )
                        if dist <= distance_threshold:
                            cluster.append(res_num)
                            unclustered.remove(res_num)
                            changed = True
                            break
            
            clusters.append(sorted(cluster))
        
        self.logger.info(f"Identified {len(clusters)} spatial clusters")
        return clusters
    
    def calculate_clustering_statistics(
        self,
        motif_residues: List[int],
        structure: Structure,
        chain_id: str = 'A',
        n_permutations: int = 1000
    ) -> Dict[str, float]:
        """
        Test if motif residues cluster spatially more than random.
        
        Uses Hopkins statistic and permutation test.
        
        Args:
            motif_residues: List of motif residue numbers
            structure: BioPython Structure object
            chain_id: PDB chain identifier
            n_permutations: Number of random permutations
            
        Returns:
            Dictionary of clustering statistics
        """
        chain = structure[0][chain_id]
        
        # Get all CA coordinates
        all_coords = {}
        for residue in chain:
            if residue.id[0] == ' ' and 'CA' in residue:
                res_num = residue.id[1]
                all_coords[res_num] = residue['CA'].get_coord()
        
        # Get motif coordinates
        motif_coords = np.array([
            all_coords[res] for res in motif_residues
            if res in all_coords
        ])
        
        if len(motif_coords) < 2:
            return {'error': 'Insufficient motif residues'}
        
        # Calculate observed mean pairwise distance
        from scipy.spatial.distance import pdist
        observed_mean_dist = np.mean(pdist(motif_coords))
        
        # Permutation test
        all_coords_array = np.array(list(all_coords.values()))
        n_motif = len(motif_coords)
        
        random_dists = []
        for _ in range(n_permutations):
            # Random sample of residues
            random_idx = np.random.choice(
                len(all_coords_array),
                size=n_motif,
                replace=False
            )
            random_coords = all_coords_array[random_idx]
            random_mean_dist = np.mean(pdist(random_coords))
            random_dists.append(random_mean_dist)
        
        random_dists = np.array(random_dists)
        
        # Calculate p-value (one-tailed: observed < random = clustered)
        p_value = np.sum(random_dists <= observed_mean_dist) / n_permutations
        
        stats = {
            'observed_mean_distance': float(observed_mean_dist),
            'expected_mean_distance': float(np.mean(random_dists)),
            'p_value': float(p_value),
            'z_score': float(
                (observed_mean_dist - np.mean(random_dists)) / np.std(random_dists)
            )
        }
        
        self.logger.info(f"Clustering statistics: {stats}")
        return stats
    
    def visualize_structure_3d(
        self,
        pdb_file: Path,
        conservation_map: Optional[Dict[int, float]] = None,
        motif_residues: Optional[List[int]] = None,
        output_html: Optional[Path] = None
    ) -> str:
        """
        Create 3D visualization of structure with conservation coloring.
        
        Args:
            pdb_file: Path to PDB file
            conservation_map: Dict mapping residue number to conservation
            motif_residues: List of motif residue numbers to highlight
            output_html: Optional path to save HTML visualization
            
        Returns:
            HTML string with embedded 3D viewer
        """
        # Read PDB file
        with open(pdb_file) as f:
            pdb_data = f.read()
        
        # Create viewer
        view = py3Dmol.view(width=800, height=600)
        view.addModel(pdb_data, 'pdb')
        
        # Default cartoon representation
        view.setStyle({'cartoon': {'color': 'lightgray'}})
        
        # Color by conservation if provided
        if conservation_map:
            for res_num, cons_score in conservation_map.items():
                # Color gradient: blue (low) -> white -> red (high)
                if cons_score < 0.5:
                    color = f'rgb({int(255*cons_score*2)},{int(255*cons_score*2)},255)'
                else:
                    color = f'rgb(255,{int(255*(2-2*cons_score))},{int(255*(2-2*cons_score))})'
                
                view.setStyle(
                    {'resi': res_num},
                    {'cartoon': {'color': color}}
                )
        
        # Highlight motif residues
        if motif_residues:
            for res_num in motif_residues:
                view.addStyle(
                    {'resi': res_num},
                    {'stick': {'color': 'red', 'radius': 0.3}}
                )
        
        view.zoomTo()
        
        # Generate HTML
        html = view._make_html()
        
        # Save if path provided
        if output_html:
            output_html.parent.mkdir(parents=True, exist_ok=True)
            with open(output_html, 'w') as f:
                f.write(html)
            self.logger.info(f"3D visualization saved to {output_html}")
        
        return html
    
    def save_structure_analysis(
        self,
        output_path: Path,
        conservation_map: Dict[int, float],
        residue_properties: Dict[int, Dict],
        motif_residues: Optional[List[int]] = None,
        clustering_stats: Optional[Dict] = None
    ):
        """
        Save comprehensive structure analysis to JSON.
        
        Args:
            output_path: Path to save JSON file
            conservation_map: Conservation scores
            residue_properties: Structural properties
            motif_residues: Motif residue numbers
            clustering_stats: Clustering statistics
        """
        data = {
            'conservation': conservation_map,
            'residue_properties': residue_properties,
            'motif_residues': motif_residues or [],
            'clustering_statistics': clustering_stats or {}
        }
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"Structure analysis saved to {output_path}")
