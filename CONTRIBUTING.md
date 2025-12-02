# Contributing to EvoMotif

First off, thank you for considering contributing to EvoMotif! It's people like you that make EvoMotif such a great tool.

## Code of Conduct

This project and everyone participating in it is governed by the [EvoMotif Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When you create a bug report, include as many details as possible:

- **Use a clear and descriptive title**
- **Describe the exact steps to reproduce the problem**
- **Provide specific examples** (code snippets, data files)
- **Describe the behavior you observed** and what you expected
- **Include error messages and stack traces**
- **Specify your environment** (OS, Python version, dependency versions)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion:

- **Use a clear and descriptive title**
- **Provide a detailed description** of the proposed functionality
- **Explain why this enhancement would be useful**
- **List any alternative solutions** you've considered

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Make your changes**:
   - Follow the [coding style](#coding-style)
   - Add tests for new functionality
   - Update documentation as needed
3. **Ensure tests pass**: `pytest tests/ -v`
4. **Ensure code is formatted**: `black evomotif tests`
5. **Create a pull request** with a clear description

## Development Setup

```bash
# Clone your fork
git clone https://github.com/yourusername/EvoMotif.git
cd EvoMotif

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install development dependencies
pip install -e ".[dev]"

# Install external tools (Ubuntu)
sudo apt-get install mafft hmmer fasttree

# Run tests
pytest tests/ -v --cov=evomotif
```

## Coding Style

We use:
- **Black** for code formatting (line length: 88)
- **flake8** for linting
- **mypy** for type checking
- **isort** for import sorting

```bash
# Format code
black evomotif tests

# Sort imports
isort evomotif tests

# Lint
flake8 evomotif

# Type check
mypy evomotif
```

### Code Style Guidelines

- Use descriptive variable names
- Add docstrings to all public functions/classes (Google style)
- Type hint function signatures
- Keep functions focused (single responsibility)
- Write tests for new code (aim for >80% coverage)

### Docstring Example

```python
def calculate_conservation(
    alignment: MultipleSeqAlignment,
    method: str = "combined"
) -> np.ndarray:
    """
    Calculate conservation scores for alignment.
    
    Args:
        alignment: Multiple sequence alignment
        method: Conservation method ('entropy', 'blosum', 'combined')
        
    Returns:
        Array of conservation scores (0-1, higher = more conserved)
        
    Raises:
        ValueError: If method is not recognized
        
    Example:
        >>> scorer = ConservationScorer()
        >>> conservation = scorer.calculate_conservation(alignment)
    """
    pass
```

## Testing

- Write unit tests for all new functions
- Use pytest fixtures for common test data
- Mock external API calls and file I/O when appropriate
- Test edge cases and error conditions

```python
# Example test
def test_motif_discovery():
    """Test basic motif discovery."""
    discoverer = MotifDiscoverer(min_conservation=0.85)
    
    # Create test data
    conservation = np.array([0.9, 0.9, 0.9, 0.5, 0.5])
    gap_freq = np.zeros(5)
    
    # Run discovery
    motifs = discoverer.scan_windows(conservation, gap_freq, window_size=3)
    
    # Assertions
    assert len(motifs) > 0
    assert motifs[0][2] >= 0.85  # mean conservation
```

## Documentation

- Update README.md for user-facing changes
- Update API documentation for new modules/functions
- Add examples for new features
- Update tutorial if workflow changes

### Building Documentation

```bash
cd docs
make html
# Open docs/_build/html/index.html
```

## Commit Messages

Follow conventional commits:

- `feat: add variant clustering analysis`
- `fix: handle empty alignment in conservation scorer`
- `docs: update installation instructions`
- `test: add tests for motif merging`
- `refactor: simplify phylogenetic inference`
- `style: format code with black`
- `chore: update dependencies`

## Versioning

We use [Semantic Versioning](https://semver.org/):
- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

## Release Process

1. Update version in `pyproject.toml` and `setup.py`
2. Update CHANGELOG.md
3. Create release tag: `git tag -a v0.2.0 -m "Release v0.2.0"`
4. Push tag: `git push origin v0.2.0`
5. GitHub Actions will build and publish to PyPI

## Questions?

- Open a [GitHub Discussion](https://github.com/yourusername/EvoMotif/discussions)
- Email: taha@example.com

## Recognition

Contributors will be:
- Listed in AUTHORS.md
- Acknowledged in release notes
- Credited in academic publications (for significant contributions)

Thank you for contributing! 🎉
