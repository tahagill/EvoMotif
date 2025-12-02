"""
Test suite for retrieval module.
"""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from evomotif.retrieval import SequenceRetriever


@pytest.fixture
def retriever():
    """Create SequenceRetriever instance."""
    return SequenceRetriever(email="test@example.com")


@pytest.fixture
def mock_sequences():
    """Create mock sequence records."""
    return [
        SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="seq1", description="Protein 1"),
        SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="seq2", description="Protein 2 fragment"),
        SeqRecord(Seq("ACDEFG"), id="seq3", description="Short protein"),
        SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWYMKL"), id="seq4", description="Protein 3"),
    ]


class TestSequenceRetriever:
    """Tests for SequenceRetriever class."""
    
    def test_initialization(self):
        """Test retriever initialization."""
        retriever = SequenceRetriever(
            email="test@example.com",
            api_key="test_key"
        )
        assert retriever is not None
    
    @patch('evomotif.retrieval.Entrez.esearch')
    @patch('evomotif.retrieval.Entrez.read')
    def test_search_protein(self, mock_read, mock_search, retriever):
        """Test protein search functionality."""
        mock_read.return_value = {"IdList": ["123", "456", "789"]}
        mock_handle = MagicMock()
        mock_search.return_value = mock_handle
        
        ids = retriever.search_protein("test_protein")
        
        assert len(ids) == 3
        assert "123" in ids
    
    def test_filter_sequences_min_length(self, retriever, mock_sequences):
        """Test sequence filtering by minimum length."""
        filtered = retriever.filter_sequences(mock_sequences, min_length=10)
        
        assert len(filtered) == 3  # seq3 is too short
        assert all(len(seq.seq) >= 10 for seq in filtered)
    
    def test_filter_sequences_fragments(self, retriever, mock_sequences):
        """Test fragment removal."""
        filtered = retriever.filter_sequences(
            mock_sequences,
            min_length=5,
            remove_fragments=True
        )
        
        # seq2 contains "fragment" in description
        assert len(filtered) == 3
        assert not any("fragment" in seq.description.lower() for seq in filtered)
    
    def test_remove_duplicates(self, retriever):
        """Test duplicate sequence removal."""
        duplicates = [
            SeqRecord(Seq("ACDEFG"), id="seq1"),
            SeqRecord(Seq("ACDEFG"), id="seq2"),  # duplicate
            SeqRecord(Seq("GHIJKL"), id="seq3"),
        ]
        
        unique = retriever.remove_duplicates(duplicates)
        
        assert len(unique) == 2
    
    def test_save_sequences(self, retriever, mock_sequences, tmp_path):
        """Test sequence saving."""
        output_path = tmp_path / "sequences.fasta"
        
        retriever.save_sequences(mock_sequences, output_path)
        
        assert output_path.exists()
        
        # Read back and verify
        from Bio import SeqIO
        loaded = list(SeqIO.parse(output_path, "fasta"))
        assert len(loaded) == len(mock_sequences)
