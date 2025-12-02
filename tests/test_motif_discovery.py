"""
Test suite for motif discovery module.
"""

import pytest
import numpy as np
from pathlib import Path
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from evomotif.motif_discovery import MotifDiscoverer, Motif


@pytest.fixture
def conservation_with_motif():
    """Create conservation array with clear motif."""
    # 100 positions, with high conservation 20-30
    conservation = np.random.uniform(0.3, 0.6, 100)
    conservation[20:30] = np.random.uniform(0.85, 0.95, 10)
    return conservation


@pytest.fixture
def gap_frequency_low():
    """Create gap frequency array with low gaps."""
    return np.random.uniform(0, 0.1, 100)


@pytest.fixture
def simple_alignment():
    """Create simple alignment."""
    sequences = [
        SeqRecord(Seq("A" * 20 + "DEFGHIKLMN" + "A" * 70), id="seq1"),
        SeqRecord(Seq("A" * 20 + "DEFGHIKLMN" + "A" * 70), id="seq2"),
        SeqRecord(Seq("A" * 20 + "DEFGHIKLMN" + "A" * 70), id="seq3"),
    ]
    return MultipleSeqAlignment(sequences)


@pytest.fixture
def discoverer():
    """Create MotifDiscoverer instance."""
    return MotifDiscoverer(
        window_sizes=[7, 9, 10],
        min_conservation=0.85,
        max_std=0.1
    )


class TestMotifDiscoverer:
    """Tests for MotifDiscoverer class."""
    
    def test_initialization(self, discoverer):
        """Test discoverer initialization."""
        assert discoverer is not None
        assert len(discoverer.window_sizes) == 3
    
    def test_scan_windows(self, discoverer, conservation_with_motif, gap_frequency_low):
        """Test sliding window scan."""
        candidates = discoverer.scan_windows(
            conservation_with_motif,
            gap_frequency_low,
            window_size=10
        )
        
        # Should find at least one candidate around position 20-30
        assert len(candidates) > 0
        
        # Check candidate format
        for start, end, mean_cons, std_cons in candidates:
            assert end - start == 10
            assert mean_cons >= discoverer.min_conservation
            assert std_cons <= discoverer.max_std
    
    def test_merge_overlapping_motifs(self, discoverer):
        """Test motif merging."""
        # Create overlapping candidates
        candidates = [
            (10, 20, 0.90, 0.05),  # Higher score
            (15, 25, 0.85, 0.08),  # Overlaps with first
            (50, 60, 0.88, 0.06),  # Separate
        ]
        
        merged = discoverer.merge_overlapping_motifs(candidates)
        
        # Should keep 2 motifs (highest scoring from overlap + separate)
        assert len(merged) == 2
        
        # Check they don't overlap
        positions = set()
        for start, end, _, _ in merged:
            motif_positions = set(range(start, end))
            assert not motif_positions & positions
            positions.update(motif_positions)
    
    def test_extract_motif_sequences(self, discoverer, simple_alignment):
        """Test motif sequence extraction."""
        sequences = discoverer.extract_motif_sequences(
            simple_alignment,
            start=20,
            end=30
        )
        
        assert len(sequences) == len(simple_alignment)
        assert all(seq == "DEFGHIKLMN" for seq in sequences)
    
    def test_discover_motifs(
        self,
        discoverer,
        simple_alignment,
        conservation_with_motif,
        gap_frequency_low
    ):
        """Test full motif discovery."""
        consensus = "A" * 20 + "DEFGHIKLMN" + "A" * 70
        
        motifs = discoverer.discover_motifs(
            simple_alignment,
            conservation_with_motif,
            gap_frequency_low,
            consensus
        )
        
        # Should find at least one motif
        assert len(motifs) > 0
        
        # Check motif properties
        for motif in motifs:
            assert isinstance(motif, Motif)
            assert motif.length > 0
            assert motif.mean_conservation >= discoverer.min_conservation
            assert len(motif.sequences) > 0
    
    def test_save_motifs_fasta(self, discoverer, tmp_path):
        """Test saving motifs to FASTA."""
        motifs = [
            Motif(
                start=10,
                end=20,
                length=10,
                mean_conservation=0.9,
                std_conservation=0.05,
                consensus_sequence="DEFGHIKLMN",
                sequences=["DEFGHIKLMN", "DEFGHIKLMN", "DEFGHIKLMN"]
            )
        ]
        
        output_dir = tmp_path / "motifs"
        discoverer.save_motifs_fasta(motifs, output_dir)
        
        # Check file was created
        motif_file = output_dir / "motif_1.fasta"
        assert motif_file.exists()
        
        # Load and verify
        from Bio import SeqIO
        records = list(SeqIO.parse(motif_file, "fasta"))
        assert len(records) == 3
