"""
Test suite for stats module.
"""

import pytest
import numpy as np
from scipy import stats

from evomotif.stats import StatisticalAnalyzer


@pytest.fixture
def analyzer():
    """Create StatisticalAnalyzer instance."""
    return StatisticalAnalyzer(random_seed=42)


class TestStatisticalAnalyzer:
    """Tests for StatisticalAnalyzer class."""
    
    def test_initialization(self, analyzer):
        """Test analyzer initialization."""
        assert analyzer is not None
    
    def test_multiple_testing_correction(self, analyzer):
        """Test multiple testing correction."""
        p_values = np.array([0.001, 0.02, 0.05, 0.1, 0.5])
        
        reject, adjusted = analyzer.multiple_testing_correction(
            p_values,
            method='fdr_bh',
            alpha=0.05
        )
        
        assert len(reject) == len(p_values)
        assert len(adjusted) == len(p_values)
        assert np.sum(reject) <= len(p_values)  # Some should be rejected
    
    def test_permutation_test(self, analyzer):
        """Test permutation test."""
        data = np.random.normal(5, 1, 100)
        observed_mean = np.mean(data)
        
        result = analyzer.permutation_test(
            observed_statistic=observed_mean,
            data=data,
            statistic_func=np.mean,
            n_permutations=1000
        )
        
        assert 'observed' in result
        assert 'expected' in result
        assert 'p_value' in result
        assert result['observed'] == observed_mean
    
    def test_bootstrap_confidence_interval(self, analyzer):
        """Test bootstrap CI calculation."""
        data = np.random.normal(10, 2, 100)
        
        lower, mean, upper = analyzer.bootstrap_confidence_interval(
            data,
            statistic_func=np.mean,
            n_bootstrap=1000,
            confidence_level=0.95
        )
        
        assert lower < mean < upper
        assert np.abs(mean - 10) < 1  # Close to true mean
    
    def test_calculate_effect_size(self, analyzer):
        """Test effect size calculation."""
        group1 = np.random.normal(10, 2, 50)
        group2 = np.random.normal(12, 2, 50)
        
        effect_sizes = analyzer.calculate_effect_size(group1, group2)
        
        assert 'cohens_d' in effect_sizes
        assert 'hedges_g' in effect_sizes
        assert 'glass_delta' in effect_sizes
        
        # Cohen's d should be negative (group1 < group2)
        assert effect_sizes['cohens_d'] < 0
    
    def test_hypergeometric_test(self, analyzer):
        """Test hypergeometric test."""
        result = analyzer.hypergeometric_test(
            overlap=15,
            motif_size=20,
            domain_size=50,
            total_size=200
        )
        
        assert 'overlap' in result
        assert 'expected' in result
        assert 'enrichment_fold' in result
        assert 'p_value' in result
        
        # With overlap > expected, should be enriched
        assert result['enrichment_fold'] > 1
    
    def test_hopkins_statistic_random(self, analyzer):
        """Test Hopkins statistic on random data."""
        # Random data should have H ~ 0.5
        random_data = np.random.uniform(0, 10, (100, 2))
        
        H = analyzer.hopkins_statistic(random_data, n_samples=20)
        
        # Should be around 0.5 for random data
        assert 0.3 < H < 0.7
    
    def test_hopkins_statistic_clustered(self, analyzer):
        """Test Hopkins statistic on clustered data."""
        # Create clustered data
        cluster1 = np.random.normal([2, 2], 0.3, (50, 2))
        cluster2 = np.random.normal([8, 8], 0.3, (50, 2))
        clustered_data = np.vstack([cluster1, cluster2])
        
        H = analyzer.hopkins_statistic(clustered_data, n_samples=20)
        
        # Should be > 0.5 for clustered data
        assert H > 0.5
    
    def test_calculate_roc_metrics(self, analyzer):
        """Test ROC metrics calculation."""
        # Create synthetic data
        np.random.seed(42)
        true_labels = np.array([0]*50 + [1]*50)
        predicted_scores = np.concatenate([
            np.random.normal(0.3, 0.2, 50),  # True negatives
            np.random.normal(0.7, 0.2, 50)   # True positives
        ])
        
        metrics = analyzer.calculate_roc_metrics(true_labels, predicted_scores)
        
        assert 'roc_auc' in metrics
        assert 'sensitivity' in metrics
        assert 'specificity' in metrics
        assert 'f1_score' in metrics
        
        # With good separation, AUC should be high
        assert metrics['roc_auc'] > 0.7
    
    def test_phylogenetic_signal_test(self, analyzer):
        """Test phylogenetic signal detection."""
        # Create distance matrix
        n = 10
        distance_matrix = np.random.uniform(0, 1, (n, n))
        distance_matrix = (distance_matrix + distance_matrix.T) / 2  # Symmetric
        np.fill_diagonal(distance_matrix, 0)
        
        # Create trait values correlated with distance
        trait_values = np.random.normal(0, 1, n)
        
        result = analyzer.phylogenetic_signal_test(
            trait_values,
            distance_matrix
        )
        
        assert 'correlation' in result
        assert 'p_value_parametric' in result
        assert 'p_value_permutation' in result
        assert 'interpretation' in result
    
    def test_test_motif_significance(self, analyzer):
        """Test motif significance testing."""
        # Create conservation array
        conservation = np.random.uniform(0.3, 0.7, 100)
        
        # Create motifs
        motifs = [
            {
                'start': 10,
                'end': 20,
                'mean_conservation': 0.9,
                'length': 10
            },
            {
                'start': 50,
                'end': 60,
                'mean_conservation': 0.85,
                'length': 10
            }
        ]
        
        results = analyzer.test_motif_significance(
            motifs,
            conservation,
            n_permutations=1000
        )
        
        assert len(results) == 2
        for result in results:
            assert 'p_value' in result
            assert 'z_score' in result
            assert 'adjusted_p_value' in result
    
    def test_save_statistical_report(self, analyzer, tmp_path):
        """Test saving statistical report."""
        output_path = tmp_path / "stats_report.json"
        
        analyzer.save_statistical_report(
            output_path,
            test1={'result': 'pass'},
            test2={'result': 'pass'}
        )
        
        assert output_path.exists()
        
        # Load and verify
        import json
        with open(output_path) as f:
            report = json.load(f)
        
        assert 'timestamp' in report
        assert 'analyses' in report
        assert 'test1' in report['analyses']
        assert 'test2' in report['analyses']
