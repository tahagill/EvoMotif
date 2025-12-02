"""
Statistical analysis module.

Provides comprehensive statistical tests and validation methods
for motif discovery, conservation analysis, and variant enrichment.
"""

import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Optional, Tuple, Callable
from scipy import stats
from scipy.stats import (
    fisher_exact, chi2_contingency, mannwhitneyu,
    ks_2samp, pearsonr, spearmanr
)
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import json
from pathlib import Path

logger = logging.getLogger(__name__)


class StatisticalAnalyzer:
    """Comprehensive statistical analysis for EvoMotif."""
    
    def __init__(self, random_seed: int = 42):
        """
        Initialize statistical analyzer.
        
        Args:
            random_seed: Random seed for reproducibility
        """
        self.logger = logging.getLogger(__name__)
        np.random.seed(random_seed)
    
    def multiple_testing_correction(
        self,
        p_values: np.ndarray,
        method: str = 'fdr_bh',
        alpha: float = 0.05
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Apply multiple testing correction.
        
        Args:
            p_values: Array of p-values
            method: Correction method ('bonferroni', 'fdr_bh', 'fdr_by')
            alpha: Significance level
            
        Returns:
            Tuple of (reject, adjusted_p_values)
        """
        reject, adjusted_p, _, _ = multipletests(
            p_values,
            alpha=alpha,
            method=method
        )
        
        self.logger.info(
            f"Multiple testing correction ({method}): "
            f"{np.sum(reject)}/{len(p_values)} significant"
        )
        
        return reject, adjusted_p
    
    def permutation_test(
        self,
        observed_statistic: float,
        data: np.ndarray,
        statistic_func: Callable,
        n_permutations: int = 10000,
        alternative: str = 'two-sided'
    ) -> Dict[str, float]:
        """
        General permutation test.
        
        Args:
            observed_statistic: Observed test statistic
            data: Data array to permute
            statistic_func: Function to calculate statistic from permuted data
            n_permutations: Number of permutations
            alternative: 'two-sided', 'greater', or 'less'
            
        Returns:
            Dictionary with test results
        """
        # Generate null distribution
        null_distribution = []
        
        for _ in range(n_permutations):
            permuted_data = np.random.permutation(data)
            null_stat = statistic_func(permuted_data)
            null_distribution.append(null_stat)
        
        null_distribution = np.array(null_distribution)
        
        # Calculate p-value
        if alternative == 'greater':
            p_value = np.sum(null_distribution >= observed_statistic) / n_permutations
        elif alternative == 'less':
            p_value = np.sum(null_distribution <= observed_statistic) / n_permutations
        else:  # two-sided
            p_value = np.sum(
                np.abs(null_distribution) >= np.abs(observed_statistic)
            ) / n_permutations
        
        # Z-score
        z_score = (observed_statistic - np.mean(null_distribution)) / np.std(null_distribution)
        
        return {
            'observed': float(observed_statistic),
            'expected': float(np.mean(null_distribution)),
            'std': float(np.std(null_distribution)),
            'p_value': float(p_value),
            'z_score': float(z_score)
        }
    
    def bootstrap_confidence_interval(
        self,
        data: np.ndarray,
        statistic_func: Callable,
        n_bootstrap: int = 10000,
        confidence_level: float = 0.95
    ) -> Tuple[float, float, float]:
        """
        Calculate bootstrap confidence interval.
        
        Args:
            data: Input data
            statistic_func: Function to calculate statistic
            n_bootstrap: Number of bootstrap samples
            confidence_level: Confidence level (e.g., 0.95 for 95% CI)
            
        Returns:
            Tuple of (lower_bound, mean, upper_bound)
        """
        bootstrap_stats = []
        
        for _ in range(n_bootstrap):
            sample = np.random.choice(data, size=len(data), replace=True)
            stat = statistic_func(sample)
            bootstrap_stats.append(stat)
        
        bootstrap_stats = np.array(bootstrap_stats)
        
        alpha = 1 - confidence_level
        lower = np.percentile(bootstrap_stats, alpha/2 * 100)
        upper = np.percentile(bootstrap_stats, (1 - alpha/2) * 100)
        mean = np.mean(bootstrap_stats)
        
        return lower, mean, upper
    
    def test_motif_significance(
        self,
        observed_motifs: List[Dict],
        alignment_conservation: np.ndarray,
        n_permutations: int = 1000
    ) -> List[Dict]:
        """
        Test if observed motifs are more conserved than random.
        
        Args:
            observed_motifs: List of motif dictionaries with 'start', 'end', 'mean_conservation'
            alignment_conservation: Full conservation array
            n_permutations: Number of random motifs to generate
            
        Returns:
            List of motifs with added significance values
        """
        results = []
        
        for motif in observed_motifs:
            start = motif['start']
            end = motif['end']
            length = end - start
            observed_cons = motif['mean_conservation']
            
            # Generate random motifs of same length
            random_cons = []
            max_start = len(alignment_conservation) - length
            
            for _ in range(n_permutations):
                random_start = np.random.randint(0, max_start)
                random_end = random_start + length
                random_mean = np.mean(alignment_conservation[random_start:random_end])
                random_cons.append(random_mean)
            
            random_cons = np.array(random_cons)
            
            # Calculate empirical p-value
            p_value = np.sum(random_cons >= observed_cons) / n_permutations
            
            # Z-score
            z_score = (observed_cons - np.mean(random_cons)) / np.std(random_cons)
            
            motif_result = motif.copy()
            motif_result.update({
                'p_value': float(p_value),
                'z_score': float(z_score),
                'expected_conservation': float(np.mean(random_cons))
            })
            
            results.append(motif_result)
        
        # Apply multiple testing correction
        p_values = np.array([m['p_value'] for m in results])
        _, adjusted_p = self.multiple_testing_correction(p_values)
        
        for i, motif_result in enumerate(results):
            motif_result['adjusted_p_value'] = float(adjusted_p[i])
        
        self.logger.info(
            f"Tested significance of {len(results)} motifs"
        )
        
        return results
    
    def calculate_effect_size(
        self,
        group1: np.ndarray,
        group2: np.ndarray
    ) -> Dict[str, float]:
        """
        Calculate effect size metrics.
        
        Args:
            group1: First group data
            group2: Second group data
            
        Returns:
            Dictionary with effect size metrics
        """
        # Cohen's d
        pooled_std = np.sqrt(
            ((len(group1) - 1) * np.var(group1, ddof=1) +
             (len(group2) - 1) * np.var(group2, ddof=1)) /
            (len(group1) + len(group2) - 2)
        )
        cohens_d = (np.mean(group1) - np.mean(group2)) / pooled_std
        
        # Hedge's g (corrected for small sample size)
        correction = 1 - (3 / (4 * (len(group1) + len(group2)) - 9))
        hedges_g = cohens_d * correction
        
        # Glass's delta (using group2 as control)
        glass_delta = (np.mean(group1) - np.mean(group2)) / np.std(group2, ddof=1)
        
        return {
            'cohens_d': float(cohens_d),
            'hedges_g': float(hedges_g),
            'glass_delta': float(glass_delta)
        }
    
    def hypergeometric_test(
        self,
        overlap: int,
        motif_size: int,
        domain_size: int,
        total_size: int
    ) -> Dict[str, float]:
        """
        Hypergeometric test for motif-domain overlap.
        
        Tests if motifs overlap with known domains more than expected by chance.
        
        Args:
            overlap: Number of overlapping positions
            motif_size: Total motif size
            domain_size: Total domain size
            total_size: Total protein length
            
        Returns:
            Dictionary with test results
        """
        # Expected overlap under null hypothesis
        expected = (motif_size * domain_size) / total_size
        
        # Hypergeometric p-value
        p_value = stats.hypergeom.sf(
            overlap - 1,  # -1 because sf is P(X > x), we want P(X >= x)
            total_size,
            domain_size,
            motif_size
        )
        
        # Enrichment fold
        enrichment = (overlap / expected) if expected > 0 else 0
        
        return {
            'overlap': overlap,
            'expected': float(expected),
            'enrichment_fold': float(enrichment),
            'p_value': float(p_value)
        }
    
    def hopkins_statistic(
        self,
        data: np.ndarray,
        n_samples: Optional[int] = None
    ) -> float:
        """
        Calculate Hopkins statistic for spatial clustering.
        
        Values close to 1 indicate strong clustering.
        Values close to 0.5 indicate random distribution.
        
        Args:
            data: N x D array of coordinates
            n_samples: Number of samples (default: 10% of data size)
            
        Returns:
            Hopkins statistic
        """
        if n_samples is None:
            n_samples = int(0.1 * len(data))
        
        n_samples = min(n_samples, len(data) - 1)
        
        # Sample real points
        sample_indices = np.random.choice(len(data), size=n_samples, replace=False)
        sampled_points = data[sample_indices]
        
        # Generate random points in same space
        mins = np.min(data, axis=0)
        maxs = np.max(data, axis=0)
        random_points = np.random.uniform(mins, maxs, size=(n_samples, data.shape[1]))
        
        # Calculate nearest neighbor distances
        from scipy.spatial.distance import cdist
        
        # Distance from sampled real points to other real points
        real_distances = cdist(sampled_points, data)
        # Exclude self-distances
        for i, idx in enumerate(sample_indices):
            real_distances[i, idx] = np.inf
        u = np.min(real_distances, axis=1)
        
        # Distance from random points to real points
        random_distances = cdist(random_points, data)
        w = np.min(random_distances, axis=1)
        
        # Hopkins statistic
        H = np.sum(w) / (np.sum(u) + np.sum(w))
        
        return float(H)
    
    def calculate_roc_metrics(
        self,
        true_labels: np.ndarray,
        predicted_scores: np.ndarray
    ) -> Dict[str, float]:
        """
        Calculate ROC curve metrics.
        
        Args:
            true_labels: Binary true labels (0/1)
            predicted_scores: Predicted scores/probabilities
            
        Returns:
            Dictionary with ROC metrics
        """
        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(true_labels, predicted_scores)
        roc_auc = auc(fpr, tpr)
        
        # Calculate precision-recall curve
        precision, recall, pr_thresholds = precision_recall_curve(
            true_labels,
            predicted_scores
        )
        pr_auc = auc(recall, precision)
        
        # Find optimal threshold (Youden's J statistic)
        j_scores = tpr - fpr
        optimal_idx = np.argmax(j_scores)
        optimal_threshold = thresholds[optimal_idx]
        
        # Calculate metrics at optimal threshold
        predictions = (predicted_scores >= optimal_threshold).astype(int)
        
        tp = np.sum((predictions == 1) & (true_labels == 1))
        tn = np.sum((predictions == 0) & (true_labels == 0))
        fp = np.sum((predictions == 1) & (true_labels == 0))
        fn = np.sum((predictions == 0) & (true_labels == 1))
        
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        f1_score = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0
        
        return {
            'roc_auc': float(roc_auc),
            'pr_auc': float(pr_auc),
            'optimal_threshold': float(optimal_threshold),
            'sensitivity': float(sensitivity),
            'specificity': float(specificity),
            'ppv': float(ppv),
            'npv': float(npv),
            'f1_score': float(f1_score)
        }
    
    def phylogenetic_signal_test(
        self,
        trait_values: np.ndarray,
        distance_matrix: np.ndarray
    ) -> Dict[str, float]:
        """
        Test for phylogenetic signal in trait data.
        
        Uses Mantel test to correlate trait similarity with phylogenetic distance.
        
        Args:
            trait_values: Trait values for each taxa
            distance_matrix: Phylogenetic distance matrix
            
        Returns:
            Dictionary with test results
        """
        # Calculate trait distance matrix
        trait_matrix = np.abs(trait_values[:, None] - trait_values[None, :])
        
        # Flatten upper triangular matrices
        n = len(trait_values)
        indices = np.triu_indices(n, k=1)
        
        trait_distances = trait_matrix[indices]
        phylo_distances = distance_matrix[indices]
        
        # Calculate correlation
        correlation, p_value = spearmanr(trait_distances, phylo_distances)
        
        # Permutation test for significance
        n_perm = 1000
        null_corr = []
        
        for _ in range(n_perm):
            perm_indices = np.random.permutation(n)
            perm_trait_matrix = trait_matrix[perm_indices][:, perm_indices]
            perm_trait_dist = perm_trait_matrix[indices]
            perm_corr, _ = spearmanr(perm_trait_dist, phylo_distances)
            null_corr.append(perm_corr)
        
        null_corr = np.array(null_corr)
        p_perm = np.sum(np.abs(null_corr) >= np.abs(correlation)) / n_perm
        
        return {
            'correlation': float(correlation),
            'p_value_parametric': float(p_value),
            'p_value_permutation': float(p_perm),
            'interpretation': 'significant phylogenetic signal' if p_perm < 0.05 else 'no phylogenetic signal'
        }
    
    def save_statistical_report(
        self,
        output_path: Path,
        **analyses: Dict
    ):
        """
        Save comprehensive statistical report.
        
        Args:
            output_path: Path to save JSON report
            **analyses: Named analysis results
        """
        report = {
            'timestamp': pd.Timestamp.now().isoformat(),
            'analyses': analyses
        }
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        self.logger.info(f"Statistical report saved to {output_path}")
