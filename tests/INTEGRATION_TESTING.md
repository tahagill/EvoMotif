# Integration Testing Guide

This directory contains comprehensive integration tests for EvoMotif that work with **real protein data** from NCBI.

## Quick Start

### 1. Basic Test (Uses Real Data)

```bash
# Test with default protein (ubiquitin)
pytest tests/test_integration.py -v -s

# Test specific protein
pytest tests/test_integration.py --protein="p53" -v -s

# Quick validation (no network/MAFFT needed)
pytest tests/test_integration.py -k quick -v
```

### 2. Command-Line Testing (Any Protein!)

```bash
# Test any protein you want
python tests/run_integration_test.py ubiquitin your.email@example.com

# With custom parameters
python tests/run_integration_test.py p53 your.email@example.com \
    --max-sequences 100 \
    --min-conservation 0.8 \
    --threads 8 \
    --output results/p53

# Quick mode (fewer sequences)
python tests/run_integration_test.py insulin your.email@example.com --quick
```

## Test Categories

### Integration Tests (`test_integration.py`)

#### 1. **Full Pipeline Tests**
Tests the complete workflow on real proteins:
- ✅ Sequence retrieval from NCBI
- ✅ Multiple sequence alignment
- ✅ Conservation scoring
- ✅ Motif discovery
- ✅ Statistical validation

**Proteins tested:**
- Ubiquitin (small, highly conserved)
- p53 (medium, cancer-related)
- Insulin (small hormone)

#### 2. **Edge Case Tests**
Tests error handling and boundary conditions:

**`test_very_short_sequences`**
- Tests with minimum length sequences
- Validates filtering logic

**`test_high_gap_alignment`**
- Tests alignment with variable-length sequences
- Validates gap handling

**`test_no_conserved_regions`**
- Tests with random sequences (should find no motifs)
- Validates false positive control

**`test_identical_sequences`**
- Tests with 100% identical sequences
- Should find entire sequence as motif

**`test_missing_data_handling`**
- Tests sequences with unknown amino acids (X)
- Validates graceful degradation

**`test_single_sequence`**
- Tests error handling with single sequence
- Should fail gracefully

**`test_duplicate_removal`**
- Tests deduplication logic
- Validates sequence uniqueness

#### 3. **Quick Validation Test**
Fast test that doesn't require network or MAFFT:
- Uses synthetic but realistic data
- Tests all core algorithms
- Perfect for CI/CD pipelines

## Usage Examples

### Test Your Own Protein

```python
# In pytest
pytest tests/test_integration.py --protein="myprotein" -v -s

# Or using the CLI
python tests/run_integration_test.py myprotein your@email.com
```

### Run Specific Test Categories

```bash
# Only integration tests
pytest tests/test_integration.py::TestIntegrationPipeline -v

# Only edge cases
pytest tests/test_integration.py::TestEdgeCases -v

# Only quick validation
pytest tests/test_integration.py::TestQuickValidation -v
```

### Debugging Failed Tests

```bash
# Run with verbose output and stop on first failure
pytest tests/test_integration.py -vsx

# Show print statements
pytest tests/test_integration.py -s

# Run specific test
pytest tests/test_integration.py::TestEdgeCases::test_no_conserved_regions -v
```

## Expected Results

### For Well-Conserved Proteins (e.g., Ubiquitin)

```
✓ Retrieved 50 sequences
✓ Alignment: 50 seqs, 76 positions, 15.2% gaps
✓ Conservation: mean=0.845, max=0.980
✓ Found 3 motifs:
  Motif 1: TLTGKTITLE (pos 5-15, cons=0.912)
  Motif 2: ESTLHLVLR (pos 25-34, cons=0.889)
  Motif 3: LIFAGKQLE (pos 45-54, cons=0.901)
✓ Validated 3 motifs:
  Motif 1: p=0.001, adj_p=0.003 **
  Motif 2: p=0.002, adj_p=0.003 **
  Motif 3: p=0.001, adj_p=0.003 **
```

### For Variable Proteins

```
✓ Conservation: mean=0.512
✓ Found 0-1 motifs (expected)
```

### For Random Sequences

```
✓ Random sequences: mean conservation=0.342, motifs=0
✓ No false positives detected
```

## Test Output Structure

```
integration_test_results/
└── protein_name/
    ├── protein_name_sequences.fasta      # Retrieved sequences
    ├── protein_name_aligned.fasta        # Aligned sequences
    ├── protein_name_conservation.json    # Conservation scores
    ├── protein_name_summary.json         # Analysis summary
    └── motifs/
        ├── motif_1.fasta                 # Motif 1 sequences
        ├── motif_2.fasta                 # Motif 2 sequences
        └── ...
```

## Requirements for Real Data Tests

### External Tools (Must be installed)
```bash
# Ubuntu/Debian
sudo apt-get install mafft hmmer fasttree

# macOS
brew install mafft hmmer fasttree
```

### Network Access
- Tests retrieve sequences from NCBI
- Requires internet connection
- Use `--quick` flag for faster testing

### Email for NCBI
- Required by NCBI Entrez API
- Use your real email address
- NCBI may block excessive requests without email

## Skipping Tests

Tests automatically skip if:
- Network unavailable (retrieval tests)
- MAFFT not installed (alignment tests)
- HMMER not installed (HMM tests)

```bash
# Force skip network tests
pytest tests/test_integration.py -m "not network"
```

## Performance

| Test Type | Sequences | Runtime | Network |
|-----------|-----------|---------|---------|
| Quick validation | 5 synthetic | ~2 sec | No |
| Edge cases | 5-20 synthetic | ~10 sec | No |
| Single protein | 20-50 real | ~5 min | Yes |
| Full suite | 100+ real | ~30 min | Yes |

## Troubleshooting

### "MAFFT not found"
```bash
# Install MAFFT
sudo apt-get install mafft  # Ubuntu
brew install mafft           # macOS
```

### "Network timeout"
```bash
# Reduce sequences or use quick mode
python tests/run_integration_test.py myprotein email@test.com --quick
```

### "No sequences found"
```bash
# Try different protein name
python tests/run_integration_test.py "p53" email@test.com

# Or increase max sequences
python tests/run_integration_test.py myprotein email@test.com --max-sequences 100
```

### "No motifs found"
This is OK! Some proteins have low conservation. Try:
- Lower conservation threshold: `--min-conservation 0.6`
- More sequences: `--max-sequences 100`
- Different protein with known motifs

## CI/CD Integration

For GitHub Actions, use quick validation:

```yaml
- name: Run integration tests
  run: |
    pytest tests/test_integration.py::TestQuickValidation -v
```

For full testing with network access:

```yaml
- name: Run full integration tests
  run: |
    pytest tests/test_integration.py --protein="ubiquitin" -v
  env:
    TEST_EMAIL: ${{ secrets.TEST_EMAIL }}
```

## Contributing

When adding new tests:
1. Add parametrized test for different proteins
2. Include edge cases
3. Add docstring explaining what's tested
4. Use `pytest.skip()` for missing dependencies
5. Log progress with print statements

## Support

If tests fail:
1. Check prerequisites (MAFFT, network, email)
2. Run with `-vsx` for detailed output
3. Check logs in test output directory
4. Open GitHub issue with error message

---

**Run your first test now:**
```bash
python tests/run_integration_test.py ubiquitin your.email@example.com --quick
```
