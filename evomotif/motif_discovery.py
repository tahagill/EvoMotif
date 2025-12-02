"""
Motif discovery module using sliding window and HMM profiling.

Identifies conserved motifs from alignment conservation scores
and builds HMM profiles for motif generalization.
"""

import numpy as np
import logging
import subprocess
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


@dataclass
class Motif:
    """Represents a discovered motif."""
    start: int  # 0-indexed start position in alignment
    end: int  # 0-indexed end position (exclusive)
    length: int
    mean_conservation: float
    std_conservation: float
    consensus_sequence: str
    sequences: List[str]  # Motif sequences from all proteins
    hmm_path: Optional[Path] = None
    evalue: Optional[float] = None


class MotifDiscoverer:
    """Discover conserved motifs using sliding window analysis."""
    
    def __init__(
        self,
        window_sizes: List[int] = [7, 9, 11, 13, 15],
        min_conservation: float = 0.85,
        max_std: float = 0.1,
        min_gap_free: float = 0.9
    ):
        """
        Initialize motif discoverer.
        
        Args:
            window_sizes: List of sliding window sizes to try
            min_conservation: Minimum mean conservation score
            max_std: Maximum std of conservation within window
            min_gap_free: Minimum fraction of sequences without gaps
        """
        self.window_sizes = window_sizes
        self.min_conservation = min_conservation
        self.max_std = max_std
        self.min_gap_free = min_gap_free
        self.logger = logging.getLogger(__name__)
    
    def scan_windows(
        self,
        conservation: np.ndarray,
        gap_frequency: np.ndarray,
        window_size: int
    ) -> List[Tuple[int, int, float, float]]:
        """
        Scan alignment with sliding window to find conserved regions.
        
        Args:
            conservation: Conservation score array
            gap_frequency: Gap frequency array
            window_size: Size of sliding window
            
        Returns:
            List of (start, end, mean_conservation, std_conservation) tuples
        """
        candidates = []
        aln_length = len(conservation)
        
        for start in range(aln_length - window_size + 1):
            end = start + window_size
            window_cons = conservation[start:end]
            window_gaps = gap_frequency[start:end]
            
            # Calculate statistics
            mean_cons = np.mean(window_cons)
            std_cons = np.std(window_cons)
            max_gap_freq = np.max(window_gaps)
            
            # Apply filters
            if (mean_cons >= self.min_conservation and
                std_cons <= self.max_std and
                max_gap_freq <= (1 - self.min_gap_free)):
                
                candidates.append((start, end, mean_cons, std_cons))
        
        return candidates
    
    def merge_overlapping_motifs(
        self,
        candidates: List[Tuple[int, int, float, float]]
    ) -> List[Tuple[int, int, float, float]]:
        """
        Merge overlapping candidate motifs.
        
        Keeps the highest-scoring motif when overlaps occur.
        
        Args:
            candidates: List of (start, end, mean_cons, std_cons)
            
        Returns:
            Merged list of non-overlapping motifs
        """
        if not candidates:
            return []
        
        # Sort by conservation score (descending)
        sorted_candidates = sorted(candidates, key=lambda x: x[2], reverse=True)
        
        merged = []
        used_positions = set()
        
        for start, end, mean_cons, std_cons in sorted_candidates:
            # Check if this motif overlaps with already selected ones
            motif_positions = set(range(start, end))
            
            if not motif_positions & used_positions:
                merged.append((start, end, mean_cons, std_cons))
                used_positions.update(motif_positions)
        
        # Sort by position
        merged.sort(key=lambda x: x[0])
        
        return merged
    
    def extract_motif_sequences(
        self,
        alignment: MultipleSeqAlignment,
        start: int,
        end: int
    ) -> List[str]:
        """
        Extract motif sequences from alignment.
        
        Args:
            alignment: MultipleSeqAlignment object
            start: Start position (0-indexed)
            end: End position (exclusive)
            
        Returns:
            List of motif sequences (without gaps)
        """
        sequences = []
        
        for record in alignment:
            motif_seq = str(record.seq[start:end]).replace('-', '')
            if motif_seq:  # Only add non-empty sequences
                sequences.append(motif_seq)
        
        return sequences
    
    def discover_motifs(
        self,
        alignment: MultipleSeqAlignment,
        conservation: np.ndarray,
        gap_frequency: np.ndarray,
        consensus: str
    ) -> List[Motif]:
        """
        Discover all motifs in alignment.
        
        Args:
            alignment: MultipleSeqAlignment object
            conservation: Conservation score array
            gap_frequency: Gap frequency array
            consensus: Consensus sequence
            
        Returns:
            List of Motif objects
        """
        all_candidates = []
        
        # Scan with each window size
        for window_size in self.window_sizes:
            candidates = self.scan_windows(
                conservation,
                gap_frequency,
                window_size
            )
            all_candidates.extend(candidates)
        
        self.logger.info(f"Found {len(all_candidates)} candidate motifs")
        
        # Merge overlapping motifs
        merged = self.merge_overlapping_motifs(all_candidates)
        
        self.logger.info(f"After merging: {len(merged)} motifs")
        
        # Create Motif objects
        motifs = []
        for idx, (start, end, mean_cons, std_cons) in enumerate(merged):
            consensus_seq = consensus[start:end].replace('-', '')
            sequences = self.extract_motif_sequences(alignment, start, end)
            
            motif = Motif(
                start=start,
                end=end,
                length=end - start,
                mean_conservation=mean_cons,
                std_conservation=std_cons,
                consensus_sequence=consensus_seq,
                sequences=sequences
            )
            
            motifs.append(motif)
            
            self.logger.info(
                f"Motif {idx+1}: pos {start+1}-{end}, "
                f"cons={mean_cons:.3f}, seq={consensus_seq}"
            )
        
        return motifs
    
    def save_motifs_fasta(
        self,
        motifs: List[Motif],
        output_dir: Path
    ):
        """
        Save motif sequences to individual FASTA files.
        
        Args:
            motifs: List of Motif objects
            output_dir: Directory to save FASTA files
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for idx, motif in enumerate(motifs):
            output_path = output_dir / f"motif_{idx+1}.fasta"
            
            records = []
            for seq_idx, seq in enumerate(motif.sequences):
                record = SeqRecord(
                    Seq(seq),
                    id=f"motif{idx+1}_seq{seq_idx+1}",
                    description=f"pos_{motif.start+1}-{motif.end}_cons_{motif.mean_conservation:.3f}"
                )
                records.append(record)
            
            SeqIO.write(records, output_path, "fasta")
            self.logger.info(f"Saved {len(records)} sequences to {output_path}")


class HMMProfileBuilder:
    """Build and search HMM profiles for motifs."""
    
    def __init__(
        self,
        hmmbuild_path: str = "hmmbuild",
        hmmsearch_path: str = "hmmsearch"
    ):
        """
        Initialize HMM profile builder.
        
        Args:
            hmmbuild_path: Path to hmmbuild executable
            hmmsearch_path: Path to hmmsearch executable
        """
        self.hmmbuild_path = hmmbuild_path
        self.hmmsearch_path = hmmsearch_path
        self.logger = logging.getLogger(__name__)
        self._check_hmmer()
    
    def _check_hmmer(self):
        """Check if HMMER tools are available."""
        for tool in [self.hmmbuild_path, self.hmmsearch_path]:
            try:
                subprocess.run(
                    [tool, "-h"],
                    capture_output=True,
                    check=True
                )
            except (FileNotFoundError, subprocess.CalledProcessError):
                raise RuntimeError(
                    f"HMMER tool not found: {tool}. "
                    "Please install HMMER: http://hmmer.org/"
                )
    
    def build_profile(
        self,
        motif_fasta: Path,
        output_hmm: Path,
        name: Optional[str] = None
    ):
        """
        Build HMM profile from motif sequences.
        
        Args:
            motif_fasta: Path to motif sequences FASTA
            output_hmm: Path to save HMM profile
            name: Optional name for the HMM
        """
        cmd = [self.hmmbuild_path]
        
        if name:
            cmd.extend(["-n", name])
        
        cmd.extend([str(output_hmm), str(motif_fasta)])
        
        self.logger.info(f"Building HMM profile: {output_hmm}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"HMM profile saved to {output_hmm}")
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"hmmbuild failed: {e.stderr}")
            raise
    
    def search_sequences(
        self,
        hmm_profile: Path,
        sequence_db: Path,
        output_table: Path,
        evalue_threshold: float = 0.01
    ) -> List[Dict]:
        """
        Search sequences using HMM profile.
        
        Args:
            hmm_profile: Path to HMM profile
            sequence_db: Path to sequence database (FASTA)
            output_table: Path to save results table
            evalue_threshold: E-value threshold for hits
            
        Returns:
            List of hit dictionaries
        """
        cmd = [
            self.hmmsearch_path,
            "--tblout", str(output_table),
            "-E", str(evalue_threshold),
            str(hmm_profile),
            str(sequence_db)
        ]
        
        self.logger.info(f"Searching sequences with HMM")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse results
            hits = self._parse_hmmsearch_table(output_table)
            self.logger.info(f"Found {len(hits)} hits (E < {evalue_threshold})")
            
            return hits
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"hmmsearch failed: {e.stderr}")
            raise
    
    def _parse_hmmsearch_table(self, table_path: Path) -> List[Dict]:
        """
        Parse hmmsearch tabular output.
        
        Args:
            table_path: Path to hmmsearch table output
            
        Returns:
            List of hit dictionaries
        """
        hits = []
        
        with open(table_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.split()
                if len(fields) < 9:
                    continue
                
                hit = {
                    'target': fields[0],
                    'query': fields[2],
                    'evalue': float(fields[4]),
                    'score': float(fields[5]),
                    'bias': float(fields[6])
                }
                hits.append(hit)
        
        return hits
    
    def build_all_profiles(
        self,
        motifs: List[Motif],
        motif_dir: Path,
        hmm_dir: Path
    ) -> List[Motif]:
        """
        Build HMM profiles for all motifs.
        
        Args:
            motifs: List of Motif objects
            motif_dir: Directory containing motif FASTA files
            hmm_dir: Directory to save HMM profiles
            
        Returns:
            Updated list of Motif objects with HMM paths
        """
        hmm_dir.mkdir(parents=True, exist_ok=True)
        
        for idx, motif in enumerate(motifs):
            motif_fasta = motif_dir / f"motif_{idx+1}.fasta"
            hmm_path = hmm_dir / f"motif_{idx+1}.hmm"
            
            if motif_fasta.exists():
                self.build_profile(
                    motif_fasta,
                    hmm_path,
                    name=f"motif_{idx+1}"
                )
                motif.hmm_path = hmm_path
        
        return motifs
