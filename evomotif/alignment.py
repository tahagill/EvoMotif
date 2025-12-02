"""
Multiple sequence alignment module.

Handles alignment using MAFFT with appropriate parameters
for motif discovery.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict, List
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

logger = logging.getLogger(__name__)


class SequenceAligner:
    """Perform multiple sequence alignment using MAFFT."""
    
    def __init__(self, mafft_path: str = "mafft"):
        """
        Initialize the aligner.
        
        Args:
            mafft_path: Path to MAFFT executable
        """
        self.mafft_path = mafft_path
        self.logger = logging.getLogger(__name__)
        self._check_mafft()
    
    def _check_mafft(self):
        """Check if MAFFT is available."""
        try:
            result = subprocess.run(
                [self.mafft_path, "--version"],
                capture_output=True,
                text=True
            )
            self.logger.info(f"MAFFT version: {result.stderr.split()[0]}")
        except FileNotFoundError:
            raise RuntimeError(
                f"MAFFT not found at {self.mafft_path}. "
                "Please install MAFFT: https://mafft.cbrc.jp/alignment/software/"
            )
    
    def align_sequences(
        self,
        input_fasta: Path,
        output_fasta: Path,
        method: str = "linsi",
        threads: int = 4,
        additional_args: Optional[List[str]] = None
    ) -> MultipleSeqAlignment:
        """
        Align sequences using MAFFT.
        
        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output alignment file
            method: MAFFT method (linsi, ginsi, einsi, auto)
            threads: Number of CPU threads
            additional_args: Additional MAFFT arguments
            
        Returns:
            MultipleSeqAlignment object
        """
        # Build command
        cmd = [self.mafft_path]
        
        # Method selection
        if method == "linsi":
            cmd.append("--localpair")
            cmd.append("--maxiterate")
            cmd.append("1000")
        elif method == "ginsi":
            cmd.append("--globalpair")
            cmd.append("--maxiterate")
            cmd.append("1000")
        elif method == "einsi":
            cmd.append("--genafpair")
            cmd.append("--maxiterate")
            cmd.append("1000")
        elif method == "auto":
            cmd.append("--auto")
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Threading
        cmd.extend(["--thread", str(threads)])
        
        # Additional arguments
        if additional_args:
            cmd.extend(additional_args)
        
        # Input file
        cmd.append(str(input_fasta))
        
        self.logger.info(f"Running MAFFT with method: {method}")
        self.logger.debug(f"Command: {' '.join(cmd)}")
        
        # Run MAFFT
        try:
            with open(output_fasta, 'w') as outfile:
                result = subprocess.run(
                    cmd,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            self.logger.info(f"Alignment saved to {output_fasta}")
            
            # Load and return alignment
            alignment = AlignIO.read(output_fasta, "fasta")
            self.logger.info(
                f"Alignment: {len(alignment)} sequences, "
                f"{alignment.get_alignment_length()} positions"
            )
            
            return alignment
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MAFFT failed: {e.stderr}")
            raise
    
    def align_with_profile(
        self,
        profile_alignment: Path,
        new_sequences: Path,
        output_fasta: Path,
        threads: int = 4
    ) -> MultipleSeqAlignment:
        """
        Align new sequences to existing profile alignment.
        
        Args:
            profile_alignment: Path to existing alignment
            new_sequences: Path to new sequences to add
            output_fasta: Path to output alignment
            threads: Number of CPU threads
            
        Returns:
            MultipleSeqAlignment object
        """
        cmd = [
            self.mafft_path,
            "--add", str(new_sequences),
            "--thread", str(threads),
            str(profile_alignment)
        ]
        
        self.logger.info("Adding sequences to profile alignment")
        
        try:
            with open(output_fasta, 'w') as outfile:
                subprocess.run(
                    cmd,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            alignment = AlignIO.read(output_fasta, "fasta")
            self.logger.info(f"Profile alignment updated: {len(alignment)} sequences")
            
            return alignment
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Profile alignment failed: {e.stderr}")
            raise
    
    def calculate_alignment_stats(
        self,
        alignment: MultipleSeqAlignment
    ) -> Dict[str, float]:
        """
        Calculate basic alignment statistics.
        
        Args:
            alignment: MultipleSeqAlignment object
            
        Returns:
            Dictionary of alignment statistics
        """
        n_seqs = len(alignment)
        aln_length = alignment.get_alignment_length()
        
        # Calculate gap statistics
        total_positions = n_seqs * aln_length
        gap_count = sum(
            str(record.seq).count('-')
            for record in alignment
        )
        gap_percentage = (gap_count / total_positions) * 100
        
        # Calculate column gap statistics
        gappy_columns = 0
        for col_idx in range(aln_length):
            column = alignment[:, col_idx]
            if column.count('-') / n_seqs > 0.5:
                gappy_columns += 1
        
        stats = {
            'n_sequences': n_seqs,
            'alignment_length': aln_length,
            'gap_percentage': gap_percentage,
            'gappy_columns': gappy_columns,
            'gappy_columns_pct': (gappy_columns / aln_length) * 100
        }
        
        self.logger.info(f"Alignment stats: {stats}")
        return stats
    
    def trim_alignment(
        self,
        alignment: MultipleSeqAlignment,
        output_fasta: Path,
        max_gap_fraction: float = 0.5,
        min_coverage: float = 0.7
    ) -> MultipleSeqAlignment:
        """
        Trim poorly aligned regions from alignment.
        
        Args:
            alignment: Input alignment
            output_fasta: Path to save trimmed alignment
            max_gap_fraction: Maximum fraction of gaps per column
            min_coverage: Minimum fraction of non-gap characters per sequence
            
        Returns:
            Trimmed MultipleSeqAlignment
        """
        aln_length = alignment.get_alignment_length()
        n_seqs = len(alignment)
        
        # Identify columns to keep
        keep_columns = []
        for col_idx in range(aln_length):
            column = alignment[:, col_idx]
            gap_fraction = column.count('-') / n_seqs
            if gap_fraction <= max_gap_fraction:
                keep_columns.append(col_idx)
        
        self.logger.info(
            f"Keeping {len(keep_columns)}/{aln_length} columns "
            f"(max gap fraction: {max_gap_fraction})"
        )
        
        # Extract kept columns
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        trimmed_records = []
        for record in alignment:
            new_seq = ''.join(str(record.seq)[i] for i in keep_columns)
            
            # Check minimum coverage
            coverage = (len(new_seq) - new_seq.count('-')) / len(new_seq)
            if coverage >= min_coverage:
                trimmed_records.append(
                    SeqRecord(Seq(new_seq), id=record.id, description="")
                )
        
        self.logger.info(
            f"Keeping {len(trimmed_records)}/{n_seqs} sequences "
            f"(min coverage: {min_coverage})"
        )
        
        # Create new alignment
        trimmed_alignment = MultipleSeqAlignment(trimmed_records)
        
        # Save
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        AlignIO.write(trimmed_alignment, output_fasta, "fasta")
        
        return trimmed_alignment
