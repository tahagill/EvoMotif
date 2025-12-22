#!/usr/bin/env python3
"""
Final integration test before PyPI release.
Tests all major features end-to-end.
"""

import sys
import tempfile
from pathlib import Path

def test_simple_api():
    """Test the simple 2-line API"""
    print("\n" + "="*70)
    print("TEST 1: Simple API (2-line usage)")
    print("="*70)
    
    import evomotif
    
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            results = evomotif.analyze_protein(
                "insulin",
                "test@email.com",
                output_dir=tmpdir,
                max_sequences=15,  # Increase for better retrieval
                verbose=False
            )
            
            print(f"✓ Analysis completed")
            print(f"✓ Sequences: {results.n_sequences}")
            print(f"✓ Motifs found: {len(results.motifs)}")
            print(f"✓ Conserved positions: {len(results.conserved_positions)}")
            print(f"✓ Output directory: {results.output_dir}")
            
            # Check files exist
            assert Path(tmpdir).exists()
            files = list(Path(tmpdir).glob("*"))
            print(f"✓ Generated {len(files)} output files")
            
            return True
        except Exception as e:
            print(f"✗ FAILED: {e}")
            return False

def test_results_object():
    """Test AnalysisResults object functionality"""
    print("\n" + "="*70)
    print("TEST 2: Results Object Features")
    print("="*70)
    
    import evomotif
    
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            results = evomotif.analyze_protein(
                "insulin",
                "test@email.com",
                output_dir=tmpdir,
                max_sequences=10,
                verbose=False
            )
            
            # Test summary
            summary = results.summary()
            assert "insulin" in summary.lower()
            print(f"✓ Summary generated ({len(summary)} chars)")
            
            # Test data access
            assert hasattr(results, 'motifs')
            assert hasattr(results, 'conserved_positions')
            assert hasattr(results, 'conservation_scores')
            print(f"✓ Data attributes accessible")
            
            # Test file export
            json_file = Path(tmpdir) / "test_export.json"
            results.export_json(json_file)
            assert json_file.exists()
            print(f"✓ JSON export working")
            
            # Test file access
            if results.files:
                for file_type, file_path in results.files.items():
                    if file_path and Path(file_path).exists():
                        print(f"  ✓ {file_type}: {Path(file_path).name}")
            
            return True
        except Exception as e:
            print(f"✗ FAILED: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_dependency_check():
    """Test dependency checking"""
    print("\n" + "="*70)
    print("TEST 3: Dependency Checking")
    print("="*70)
    
    try:
        from evomotif import EvoMotifPipeline
        
        pipeline = EvoMotifPipeline()
        deps = pipeline.check_dependencies()
        
        for tool, available in deps.items():
            status = "✓" if available else "✗"
            print(f"{status} {tool}: {'available' if available else 'NOT FOUND'}")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        return False

def test_imports():
    """Test all main imports work"""
    print("\n" + "="*70)
    print("TEST 4: Module Imports")
    print("="*70)
    
    try:
        # Test main imports
        import evomotif
        print("✓ evomotif")
        
        from evomotif import analyze_protein, EvoMotifPipeline, AnalysisResults
        print("✓ High-level API (analyze_protein, EvoMotifPipeline, AnalysisResults)")
        
        # Test individual modules
        from evomotif.retrieval import SequenceRetriever
        print("✓ retrieval.SequenceRetriever")
        
        from evomotif.alignment import SequenceAligner
        print("✓ alignment.SequenceAligner")
        
        from evomotif.conservation import ConservationScorer
        print("✓ conservation.ConservationScorer")
        
        from evomotif.motif_discovery import MotifDiscoverer
        print("✓ motif_discovery.MotifDiscoverer")
        
        from evomotif.stats import StatisticalAnalyzer
        print("✓ stats.StatisticalAnalyzer")
        
        from evomotif.phylogeny import PhylogeneticAnalyzer
        print("✓ phylogeny.PhylogeneticAnalyzer")
        
        from evomotif.structure import StructureMapper
        print("✓ structure.StructureMapper")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_module_attributes():
    """Test module has correct attributes"""
    print("\n" + "="*70)
    print("TEST 5: Module Attributes")
    print("="*70)
    
    try:
        import evomotif
        
        # Check __all__
        assert hasattr(evomotif, '__all__')
        print(f"✓ __all__ defined ({len(evomotif.__all__)} exports)")
        
        # Check version
        assert hasattr(evomotif, '__version__')
        print(f"✓ __version__: {evomotif.__version__}")
        
        # Check main functions are exported
        assert 'analyze_protein' in evomotif.__all__
        assert 'EvoMotifPipeline' in evomotif.__all__
        assert 'AnalysisResults' in evomotif.__all__
        print(f"✓ Main API functions exported")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        return False

def test_conservation_methods():
    """Test conservation calculation methods"""
    print("\n" + "="*70)
    print("TEST 6: Conservation Methods")
    print("="*70)
    
    try:
        from evomotif.conservation import ConservationScorer
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        scorer = ConservationScorer()
        
        # Create simple alignment
        seqs = [
            SeqRecord(Seq("MAKKDKLVVGVS"), id="seq1"),
            SeqRecord(Seq("MAKKDKLVVGVS"), id="seq2"),
            SeqRecord(Seq("MAKKDKLVVGVS"), id="seq3"),
        ]
        alignment = MultipleSeqAlignment(seqs)
        
        # Test Shannon entropy
        shannon = scorer.calculate_shannon_entropy(alignment)
        print(f"✓ Shannon entropy: {len(shannon)} positions")
        
        # Test BLOSUM score
        blosum = scorer.calculate_blosum_score(alignment)
        print(f"✓ BLOSUM62 score: {len(blosum)} positions")
        
        # Test combined
        combined = scorer.calculate_conservation_scores(alignment)
        print(f"✓ Combined conservation: {len(combined)} positions")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_alphafold_feature():
    """Test AlphaFold confidence extraction"""
    print("\n" + "="*70)
    print("TEST 7: AlphaFold Integration")
    print("="*70)
    
    try:
        from evomotif.structure import StructureMapper
        from Bio.PDB import Structure, Model, Chain, Residue, Atom
        
        mapper = StructureMapper()
        
        # Create mock AlphaFold structure
        structure = Structure.Structure('test')
        model = Model.Model(0)
        chain = Chain.Chain('A')
        
        for i in range(1, 11):
            residue = Residue.Residue((' ', i, ' '), 'ALA', i)
            atom = Atom.Atom('CA', [0, 0, 0], 50.0 + i*5, 1.0, ' ', 'CA', i, 'C')
            residue.add(atom)
            chain.add(residue)
        
        model.add(chain)
        structure.add(model)
        
        # Extract confidence
        confidence = mapper.get_alphafold_confidence(structure, chain_id='A')
        
        assert len(confidence) == 10
        print(f"✓ AlphaFold confidence extraction: {len(confidence)} residues")
        print(f"  Sample pLDDT scores: {list(confidence.values())[:5]}")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("\n" + "#"*70)
    print("# EVOMOTIF FINAL INTEGRATION TEST")
    print("# Testing all features before PyPI release")
    print("#"*70)
    
    tests = [
        ("Imports", test_imports),
        ("Module Attributes", test_module_attributes),
        ("Dependency Check", test_dependency_check),
        ("Conservation Methods", test_conservation_methods),
        ("AlphaFold Integration", test_alphafold_feature),
        ("Simple API", test_simple_api),
        ("Results Object", test_results_object),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"\n✗ {name} CRASHED: {e}")
            results.append((name, False))
    
    # Summary
    print("\n" + "#"*70)
    print("# TEST SUMMARY")
    print("#"*70)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{status}: {name}")
    
    print("\n" + "="*70)
    print(f"RESULTS: {passed}/{total} tests passed")
    print("="*70)
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! Ready for PyPI release! 🚀")
        return 0
    else:
        print(f"\n❌ {total - passed} test(s) failed. Fix issues before release.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
