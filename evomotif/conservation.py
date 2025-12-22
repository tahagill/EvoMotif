"""
Conservation scoring module.

Implements Shannon entropy and BLOSUM62-weighted conservation scoring
for multiple sequence alignments.
"""

import numpy as np
import logging
from typing import Dict, List, Tuple, Optional
from Bio.Align import MultipleSeqAlignment
from Bio.Align import substitution_matrices
from pathlib import Path
import json

logger = logging.getLogger(__name__)


class ConservationScorer:
    """Calculate conservation scores for protein alignments."""
    
    # Standard amino acids
    AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
    
    def __init__(self):
        """Initialize conservation scorer."""
        self.logger = logging.getLogger(__name__)
        self.blosum62 = substitution_matrices.load("BLOSUM62")
    
    def calculate_shannon_entropy(
        self,
        alignment: MultipleSeqAlignment,
        ignore_gaps: bool = True,
        pseudocount: float = 0.0001
    ) -> np.ndarray:
        """
        Calculate Shannon entropy for each alignment column.
        
        Lower entropy = higher conservation
        
        Args:
            alignment: MultipleSeqAlignment object
            ignore_gaps: Exclude gaps from entropy calculation
            pseudocount: Small value to avoid log(0)
            
        Returns:
            Array of entropy values (one per column)
        """
        n_seqs = len(alignment)
        aln_length = alignment.get_alignment_length()
        entropies = np.zeros(aln_length)
        
        for col_idx in range(aln_length):
            column = alignment[:, col_idx].upper()
            
            # Count amino acid frequencies
            aa_counts = {}
            total = 0
            
            for aa in column:
                if ignore_gaps and aa == '-':
                    continue
                if aa in self.AMINO_ACIDS:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1
                    total += 1
            
            if total == 0:
                entropies[col_idx] = np.nan
                continue
            
            # Calculate Shannon entropy
            entropy = 0.0
            for aa, count in aa_counts.items():
                p = (count + pseudocount) / (total + pseudocount * 20)
                entropy -= p * np.log2(p)
            
            entropies[col_idx] = entropy
        
        return entropies
    
    def calculate_blosum_score(
        self,
        alignment: MultipleSeqAlignment,
        ignore_gaps: bool = True
    ) -> np.ndarray:
        """
        Calculate BLOSUM62-weighted conservation score.
        
        For each column, calculates average BLOSUM62 score
        of all amino acids against the consensus.
        
        Args:
            alignment: MultipleSeqAlignment object
            ignore_gaps: Exclude gaps from calculation
            
        Returns:
            Array of BLOSUM scores (one per column)
        """
        aln_length = alignment.get_alignment_length()
        blosum_scores = np.zeros(aln_length)
        
        for col_idx in range(aln_length):
            column = alignment[:, col_idx].upper()
            
            # Get amino acids (excluding gaps if specified)
            aas = [aa for aa in column if aa in self.AMINO_ACIDS]
            
            if not aas:
                blosum_scores[col_idx] = np.nan
                continue
            
            # Find consensus (most frequent amino acid)
            from collections import Counter
            consensus = Counter(aas).most_common(1)[0][0]
            
            # Calculate average BLOSUM score against consensus
            scores = []
            for aa in aas:
                try:
                    score = self.blosum62[aa, consensus]
                    scores.append(score)
                except KeyError:
                    continue
            
            if scores:
                blosum_scores[col_idx] = np.mean(scores)
            else:
                blosum_scores[col_idx] = np.nan
        
        return blosum_scores
    
    def calculate_combined_conservation(
        self,
        alignment: MultipleSeqAlignment,
        entropy_weight: float = 0.5,
        blosum_weight: float = 0.5
    ) -> np.ndarray:
        """
        Calculate combined conservation score.
        
        Combines normalized Shannon entropy and BLOSUM scores.
        Higher values = higher conservation.
        
        Args:
            alignment: MultipleSeqAlignment object
            entropy_weight: Weight for entropy component
            blosum_weight: Weight for BLOSUM component
            
        Returns:
            Array of conservation scores in [0, 1]
        """
        # Calculate components
        entropy = self.calculate_shannon_entropy(alignment)
        blosum = self.calculate_blosum_score(alignment)
        
        # Normalize entropy (invert so high = conserved)
        # Max entropy for 20 amino acids = log2(20) ≈ 4.32
        max_entropy = np.log2(20)
        entropy_normalized = 1 - (entropy / max_entropy)
        entropy_normalized = np.clip(entropy_normalized, 0, 1)
        
        # Normalize BLOSUM (shift and scale to [0, 1])
        # BLOSUM62 range is approximately [-4, 11]
        blosum_normalized = (blosum + 4) / 15
        blosum_normalized = np.clip(blosum_normalized, 0, 1)
        
        # Combine
        conservation = (
            entropy_weight * entropy_normalized +
            blosum_weight * blosum_normalized
        )
        
        # Handle NaN values (gap-only columns)
        conservation = np.nan_to_num(conservation, nan=0.0)
        
        return conservation
    
    def calculate_conservation_scores(
        self,
        alignment: MultipleSeqAlignment,
        method: str = "combined",
        weights: tuple = (0.5, 0.5)
    ) -> np.ndarray:
        """
        Calculate conservation scores using specified method.
        
        Args:
            alignment: MultipleSeqAlignment object
            method: 'shannon', 'blosum62', or 'combined'
            weights: (entropy_weight, blosum_weight) if method='combined'
            
        Returns:
            Array of conservation scores
        """
        if method == "shannon":
            entropy = self.calculate_shannon_entropy(alignment)
            max_entropy = np.log2(20)
            return 1 - (entropy / max_entropy)
        elif method == "blosum62":
            blosum = self.calculate_blosum_score(alignment)
            return (blosum + 4) / 15
        else:  # combined
            return self.calculate_combined_conservation(
                alignment,
                entropy_weight=weights[0],
                blosum_weight=weights[1]
            )
    
    def calculate_gap_frequency(
        self,
        alignment: MultipleSeqAlignment
    ) -> np.ndarray:
        """
        Calculate fraction of gaps in each column.
        
        Args:
            alignment: MultipleSeqAlignment object
            
        Returns:
            Array of gap frequencies [0, 1]
        """
        n_seqs = len(alignment)
        aln_length = alignment.get_alignment_length()
        gap_freq = np.zeros(aln_length)
        
        for col_idx in range(aln_length):
            column = alignment[:, col_idx]
            gap_freq[col_idx] = column.count('-') / n_seqs
        
        return gap_freq
    
    def calculate_sequence_identity(
        self,
        alignment: MultipleSeqAlignment
    ) -> np.ndarray:
        """
        Calculate percentage identity for each column.
        
        Identity = fraction of sequences matching the consensus amino acid.
        Gaps are not counted as matching the consensus (treated as non-identity).
        
        Args:
            alignment: MultipleSeqAlignment object
            
        Returns:
            Array of identity percentages [0, 1]
        """
        n_seqs = len(alignment)
        aln_length = alignment.get_alignment_length()
        identity = np.zeros(aln_length)
        
        for col_idx in range(aln_length):
            column = alignment[:, col_idx].upper()
            
            # Get most common amino acid
            from collections import Counter
            counts = Counter(aa for aa in column if aa in self.AMINO_ACIDS)
            
            if counts:
                max_count = counts.most_common(1)[0][1]
                identity[col_idx] = max_count / n_seqs
            else:
                identity[col_idx] = 0.0
        
        return identity
    
    def get_consensus_sequence(
        self,
        alignment: MultipleSeqAlignment,
        threshold: float = 0.5
    ) -> str:
        """
        Generate consensus sequence from alignment.
        
        Args:
            alignment: MultipleSeqAlignment object
            threshold: Minimum frequency for consensus (else 'X')
            
        Returns:
            Consensus sequence string
        """
        aln_length = alignment.get_alignment_length()
        consensus = []
        
        for col_idx in range(aln_length):
            column = alignment[:, col_idx].upper()
            
            # Count amino acids
            from collections import Counter
            counts = Counter(aa for aa in column if aa in self.AMINO_ACIDS)
            
            if counts:
                most_common, count = counts.most_common(1)[0]
                freq = count / len(column)
                
                if freq >= threshold:
                    consensus.append(most_common)
                else:
                    consensus.append('X')
            else:
                consensus.append('-')
        
        return ''.join(consensus)
    
    def save_conservation_scores(
        self,
        alignment: MultipleSeqAlignment,
        output_path: Path,
        include_stats: bool = True
    ):
        """
        Calculate and save all conservation metrics.
        
        Args:
            alignment: MultipleSeqAlignment object
            output_path: Path to save JSON file
            include_stats: Include summary statistics
        """
        # Calculate all metrics
        entropy = self.calculate_shannon_entropy(alignment)
        blosum = self.calculate_blosum_score(alignment)
        conservation = self.calculate_combined_conservation(alignment)
        gap_freq = self.calculate_gap_frequency(alignment)
        identity = self.calculate_sequence_identity(alignment)
        consensus = self.get_consensus_sequence(alignment)
        
        # Prepare data
        data = {
            'alignment_length': alignment.get_alignment_length(),
            'n_sequences': len(alignment),
            'consensus': consensus,
            'positions': []
        }
        
        # Add per-position data
        for i in range(alignment.get_alignment_length()):
            pos_data = {
                'position': i + 1,
                'consensus': consensus[i],
                'entropy': float(entropy[i]) if not np.isnan(entropy[i]) else None,
                'blosum_score': float(blosum[i]) if not np.isnan(blosum[i]) else None,
                'conservation': float(conservation[i]),
                'gap_frequency': float(gap_freq[i]),
                'identity': float(identity[i])
            }
            data['positions'].append(pos_data)
        
        # Add summary statistics
        if include_stats:
            data['statistics'] = {
                'mean_conservation': float(np.nanmean(conservation)),
                'median_conservation': float(np.nanmedian(conservation)),
                'std_conservation': float(np.nanstd(conservation)),
                'mean_entropy': float(np.nanmean(entropy)),
                'mean_gap_frequency': float(np.mean(gap_freq)),
                'mean_identity': float(np.mean(identity))
            }
        
        # Save
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"Conservation scores saved to {output_path}")
