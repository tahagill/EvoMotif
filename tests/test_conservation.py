"""
Test suite for conservation module.
"""

import pytest
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from evomotif.conservation import ConservationScorer


@pytest.fixture
def simple_alignment():
    """Create a simple test alignment."""
    sequences = [
        SeqRecord(Seq("ACDEFG"), id="seq1"),
        SeqRecord(Seq("ACDEFG"), id="seq2"),
        SeqRecord(Seq("ACDEFG"), id="seq3"),
        SeqRecord(Seq("ACXEFG"), id="seq4"),  # One mutation
    ]
    return MultipleSeqAlignment(sequences)


@pytest.fixture
def variable_alignment():
    """Create alignment with variable conservation."""
    sequences = [
        SeqRecord(Seq("AAAAAVVVVV"), id="seq1"),
        SeqRecord(Seq("AAAAAVVVVV"), id="seq2"),
        SeqRecord(Seq("AAAAAXXXXX"), id="seq3"),
        SeqRecord(Seq("AAAAAYYYYY"), id="seq4"),
    ]
    return MultipleSeqAlignment(sequences)


@pytest.fixture
def scorer():
    """Create ConservationScorer instance."""
    return ConservationScorer()


class TestConservationScorer:
    """Tests for ConservationScorer class."""
    
    def test_initialization(self, scorer):
        """Test scorer initialization."""
        assert scorer is not None
        assert scorer.blosum62 is not None
    
    def test_shannon_entropy_perfect_conservation(self, simple_alignment, scorer):
        """Test Shannon entropy on perfectly conserved columns."""
        entropy = scorer.calculate_shannon_entropy(simple_alignment)
        
        # First 3 positions are identical (low entropy)
        assert entropy[0] < 0.5
        assert entropy[1] < 0.5
        assert entropy[2] < 0.5
    
    def test_shannon_entropy_variable_position(self, variable_alignment, scorer):
        """Test Shannon entropy on variable columns."""
        entropy = scorer.calculate_shannon_entropy(variable_alignment)
        
        # First 5 positions are conserved
        assert np.mean(entropy[:5]) < 0.5
        
        # Last 5 positions are variable
        assert np.mean(entropy[5:]) > 0.5
    
    def test_blosum_score(self, simple_alignment, scorer):
        """Test BLOSUM scoring."""
        blosum = scorer.calculate_blosum_score(simple_alignment)
        
        assert len(blosum) == simple_alignment.get_alignment_length()
        assert not np.all(np.isnan(blosum))
    
    def test_combined_conservation(self, simple_alignment, scorer):
        """Test combined conservation score."""
        conservation = scorer.calculate_combined_conservation(simple_alignment)
        
        assert len(conservation) == simple_alignment.get_alignment_length()
        assert np.all((conservation >= 0) & (conservation <= 1))
    
    def test_gap_frequency(self, scorer):
        """Test gap frequency calculation."""
        sequences = [
            SeqRecord(Seq("A-C"), id="seq1"),
            SeqRecord(Seq("A-C"), id="seq2"),
            SeqRecord(Seq("AAC"), id="seq3"),
        ]
        alignment = MultipleSeqAlignment(sequences)
        
        gap_freq = scorer.calculate_gap_frequency(alignment)
        
        assert gap_freq[0] == 0.0  # No gaps
        assert gap_freq[1] == 2/3  # 2 out of 3 have gaps
        assert gap_freq[2] == 0.0  # No gaps
    
    def test_sequence_identity(self, simple_alignment, scorer):
        """Test sequence identity calculation."""
        identity = scorer.calculate_sequence_identity(simple_alignment)
        
        # Most positions are 100% identical
        assert np.mean(identity) > 0.9
    
    def test_consensus_sequence(self, simple_alignment, scorer):
        """Test consensus sequence generation."""
        consensus = scorer.get_consensus_sequence(simple_alignment, threshold=0.5)
        
        assert len(consensus) == simple_alignment.get_alignment_length()
        assert consensus == "ACDEFG"
    
    def test_save_conservation_scores(self, simple_alignment, scorer, tmp_path):
        """Test saving conservation scores."""
        output_path = tmp_path / "conservation.json"
        
        scorer.save_conservation_scores(simple_alignment, output_path)
        
        assert output_path.exists()
        
        # Load and verify
        import json
        with open(output_path) as f:
            data = json.load(f)
        
        assert 'alignment_length' in data
        assert 'n_sequences' in data
        assert 'positions' in data
        assert len(data['positions']) == simple_alignment.get_alignment_length()
