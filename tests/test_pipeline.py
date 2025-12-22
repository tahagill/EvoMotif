"""
Tests for high-level pipeline API.
"""

import pytest
from pathlib import Path
import json

from evomotif import analyze_protein, EvoMotifPipeline, AnalysisResults


class TestPipelineAPI:
    """Test high-level pipeline interface."""
    
    def test_pipeline_initialization(self):
        """Test pipeline can be initialized."""
        pipeline = EvoMotifPipeline()
        assert pipeline is not None
    
    def test_check_dependencies(self):
        """Test dependency checking."""
        pipeline = EvoMotifPipeline()
        deps = pipeline.check_dependencies(['mafft', 'fasttree'])
        
        assert isinstance(deps, dict)
        assert 'mafft' in deps
        assert 'fasttree' in deps
        assert isinstance(deps['mafft'], bool)
    
    def test_analyze_protein_basic(self, tmp_path):
        """Test basic protein analysis."""
        # Skip if no internet or NCBI access
        pytest.skip("Requires NCBI access - run manually")
        
        results = analyze_protein(
            protein_name="ubiquitin",
            email="test@example.com",
            output_dir=str(tmp_path / "test_output"),
            max_sequences=10,
            check_deps=False  # Skip dependency check in tests
        )
        
        assert isinstance(results, AnalysisResults)
        assert results.protein == "ubiquitin"
        assert results.n_sequences > 0
        assert results.output_dir.exists()
    
    def test_invalid_parameters(self):
        """Test that invalid parameters raise errors."""
        with pytest.raises(ValueError):
            analyze_protein(
                protein_name="",
                email="test@example.com"
            )
        
        with pytest.raises(ValueError):
            analyze_protein(
                protein_name="test",
                email=""
            )
        
        with pytest.raises(ValueError):
            analyze_protein(
                protein_name="test",
                email="test@example.com",
                min_conservation=1.5  # Invalid range
            )
        
        with pytest.raises(ValueError):
            analyze_protein(
                protein_name="test",
                email="test@example.com",
                max_sequences=2  # Too few
            )
    
    def test_analysis_results_object(self):
        """Test AnalysisResults container."""
        test_data = {
            'protein': 'test_protein',
            'n_sequences': 10,
            'motifs': [
                {'start': 0, 'end': 5, 'sequence': 'ACDEF', 'conservation': 0.8}
            ],
            'conserved_positions': [
                {'position': 0, 'amino_acid': 'A', 'conservation': 0.9}
            ],
            'conservation': [0.8, 0.9, 0.7],
            'output_dir': '/tmp/test',
            'files': {'sequences': '/tmp/test/sequences.fasta'}
        }
        
        results = AnalysisResults(test_data)
        
        assert results.protein == 'test_protein'
        assert results.n_sequences == 10
        assert len(results.motifs) == 1
        assert len(results.conserved_positions) == 1
        
        # Test summary generation
        summary = results.summary()
        assert 'test_protein' in summary
        assert '10' in summary
        
        # Test repr
        repr_str = repr(results)
        assert 'test_protein' in repr_str
    
    def test_results_export_json(self, tmp_path):
        """Test JSON export."""
        test_data = {
            'protein': 'test',
            'n_sequences': 5,
            'motifs': [],
            'conserved_positions': [],
            'output_dir': str(tmp_path),
            'files': {}
        }
        
        results = AnalysisResults(test_data)
        
        json_file = tmp_path / "export.json"
        results.export_json(json_file)
        
        assert json_file.exists()
        
        # Verify content
        with open(json_file) as f:
            data = json.load(f)
        
        assert data['protein'] == 'test'
        assert data['n_sequences'] == 5
    
    def test_results_get_file(self):
        """Test file path retrieval."""
        test_data = {
            'protein': 'test',
            'n_sequences': 5,
            'motifs': [],
            'conserved_positions': [],
            'output_dir': '/tmp/test',
            'files': {
                'sequences': Path('/tmp/test/sequences.fasta'),
                'alignment': Path('/tmp/test/alignment.fasta')
            }
        }
        
        results = AnalysisResults(test_data)
        
        seq_file = results.get_file('sequences')
        assert seq_file == Path('/tmp/test/sequences.fasta')
        
        missing = results.get_file('nonexistent')
        assert missing is None


class TestIntegrationWithPipeline:
    """Integration tests using the pipeline API."""
    
    def test_pipeline_with_mock_data(self, tmp_path):
        """Test pipeline with minimal mock data."""
        # This would require mocking NCBI, MAFFT, etc.
        # For now, just verify the structure
        pytest.skip("Requires full environment setup")
    
    def test_convenience_function(self):
        """Test that convenience function exists and is callable."""
        assert callable(analyze_protein)
        
        # Test that it requires correct parameters
        with pytest.raises((ValueError, TypeError)):
            analyze_protein()  # No parameters


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
