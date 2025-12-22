"""
Tests for structure mapping module.
"""

import pytest
import numpy as np
from pathlib import Path
from Bio.PDB import Structure, Model, Chain, Residue
from Bio.PDB.Atom import Atom

from evomotif.structure import StructureMapper


def create_test_structure(plddt_scores=None):
    """
    Create a test structure with optional pLDDT scores in B-factors.
    
    Args:
        plddt_scores: List of (res_num, plddt) tuples, or None for regular structure
    """
    structure = Structure.Structure('test')
    model = Model.Model(0)
    chain = Chain.Chain('A')
    
    if plddt_scores is None:
        # Regular structure with default B-factors
        plddt_scores = [(1, 50.0), (2, 50.0), (3, 50.0)]
    
    for res_num, bfactor in plddt_scores:
        residue = Residue.Residue((' ', res_num, ' '), 'ALA', '')
        ca_atom = Atom(
            name='CA',
            coord=np.array([0.0, 0.0, 0.0]),
            bfactor=bfactor,
            occupancy=1.0,
            altloc=' ',
            fullname=' CA ',
            serial_number=res_num,
            element='C'
        )
        residue.add(ca_atom)
        chain.add(residue)
    
    model.add(chain)
    structure.add(model)
    return structure


class TestStructureMapper:
    """Test structure mapping functionality."""
    
    def test_initialization(self):
        """Test mapper initialization."""
        mapper = StructureMapper()
        assert mapper is not None
        assert mapper.pdb_file_path is None
    
    def test_alphafold_confidence_extraction(self):
        """Test extraction of pLDDT scores from AlphaFold models."""
        # Create mock AlphaFold structure with pLDDT scores
        plddt_data = [
            (1, 95.5),   # Very high confidence
            (2, 85.2),   # High confidence
            (3, 72.8),   # Confident
            (4, 58.3),   # Low confidence
            (5, 42.1),   # Very low confidence
        ]
        structure = create_test_structure(plddt_data)
        
        # Extract confidence
        mapper = StructureMapper()
        confidence = mapper.get_alphafold_confidence(structure, chain_id='A')
        
        # Verify all residues present
        assert len(confidence) == 5
        
        # Verify values
        for res_num, expected_plddt in plddt_data:
            assert res_num in confidence
            assert abs(confidence[res_num] - expected_plddt) < 0.01
    
    def test_alphafold_confidence_empty_chain(self):
        """Test confidence extraction with empty chain."""
        structure = Structure.Structure('test')
        model = Model.Model(0)
        chain = Chain.Chain('A')
        model.add(chain)
        structure.add(model)
        
        mapper = StructureMapper()
        confidence = mapper.get_alphafold_confidence(structure, chain_id='A')
        
        assert confidence == {}
    
    def test_alphafold_confidence_missing_ca(self):
        """Test confidence extraction when CA atoms are missing."""
        structure = Structure.Structure('test')
        model = Model.Model(0)
        chain = Chain.Chain('A')
        
        # Residue without CA atom
        residue = Residue.Residue((' ', 1, ' '), 'ALA', '')
        cb_atom = Atom(
            name='CB',
            coord=np.array([0.0, 0.0, 0.0]),
            bfactor=95.0,
            occupancy=1.0,
            altloc=' ',
            fullname=' CB ',
            serial_number=1,
            element='C'
        )
        residue.add(cb_atom)
        chain.add(residue)
        
        model.add(chain)
        structure.add(model)
        
        mapper = StructureMapper()
        confidence = mapper.get_alphafold_confidence(structure, chain_id='A')
        
        # Should return empty since no CA atoms
        assert confidence == {}
    
    def test_alphafold_confidence_ranges(self):
        """Test that extracted pLDDT scores cover expected ranges."""
        plddt_data = [
            (1, 98.5),   # Very high (>90)
            (2, 75.0),   # Confident (70-90)
            (3, 55.0),   # Low (50-70)
            (4, 30.0),   # Very low (<50)
        ]
        structure = create_test_structure(plddt_data)
        
        mapper = StructureMapper()
        confidence = mapper.get_alphafold_confidence(structure, chain_id='A')
        
        # Check ranges
        assert confidence[1] > 90  # Very high
        assert 70 <= confidence[2] < 90  # Confident
        assert 50 <= confidence[3] < 70  # Low
        assert confidence[4] < 50  # Very low
    
    def test_map_conservation_to_structure(self):
        """Test mapping conservation scores to structure residues."""
        mapper = StructureMapper()
        
        # Conservation scores for alignment positions
        conservation = np.array([0.8, 0.9, 0.7, 0.6])
        
        # Mapping from alignment to structure
        alignment_to_structure = {
            0: 10,  # Alignment pos 0 -> Structure res 10
            1: 11,  # Alignment pos 1 -> Structure res 11
            2: 12,  # Alignment pos 2 -> Structure res 12
            3: 13,  # Alignment pos 3 -> Structure res 13
        }
        
        # Map conservation
        conservation_map = mapper.map_conservation_to_structure(
            conservation,
            alignment_to_structure
        )
        
        # Verify mapping
        assert len(conservation_map) == 4
        assert conservation_map[10] == 0.8
        assert conservation_map[11] == 0.9
        assert conservation_map[12] == 0.7
        assert conservation_map[13] == 0.6
    
    def test_map_conservation_partial_overlap(self):
        """Test conservation mapping with partial overlap."""
        mapper = StructureMapper()
        
        conservation = np.array([0.8, 0.9, 0.7])
        
        # Some alignment positions don't map to structure
        alignment_to_structure = {
            0: 10,
            2: 12,  # Position 1 missing
        }
        
        conservation_map = mapper.map_conservation_to_structure(
            conservation,
            alignment_to_structure
        )
        
        # Only mapped positions should be present
        assert len(conservation_map) == 2
        assert conservation_map[10] == 0.8
        assert conservation_map[12] == 0.7
        assert 11 not in conservation_map


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
