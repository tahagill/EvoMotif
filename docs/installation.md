# Installation Guide

## Requirements

### Python Requirements

- Python 3.8 or higher
- pip or conda package manager

### System Requirements

- **RAM**: 8GB minimum, 16GB recommended
- **Disk Space**: 10GB minimum for software and data
- **OS**: Linux (recommended), macOS, or Windows with WSL

## Installation Methods

### Method 1: pip install (Recommended)

```bash
# Install from PyPI (when published)
pip install evomotif

# Or install from GitHub
pip install git+https://github.com/yourusername/EvoMotif.git
```

### Method 2: Development Installation

```bash
# Clone repository
git clone https://github.com/yourusername/EvoMotif.git
cd EvoMotif

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"
```

### Method 3: Conda Installation

```bash
# Create conda environment
conda create -n evomotif python=3.10
conda activate evomotif

# Install from conda-forge (when available)
conda install -c conda-forge evomotif

# Or install with pip in conda environment
pip install evomotif
```

## External Dependencies

EvoMotif requires several external bioinformatics tools.

### MAFFT (Required)

Multiple sequence alignment tool.

**Ubuntu/Debian:**
```bash
sudo apt-get install mafft
```

**macOS:**
```bash
brew install mafft
```

**From source:**
```bash
wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz
tar xfz mafft-7.505-with-extensions-src.tgz
cd mafft-7.505-with-extensions/core/
make clean
make
sudo make install
```

### HMMER (Required)

Hidden Markov Model tools for profile HMMs.

**Ubuntu/Debian:**
```bash
sudo apt-get install hmmer
```

**macOS:**
```bash
brew install hmmer
```

**From source:**
```bash
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.3.2
./configure
make
sudo make install
```

### FastTree (Required for phylogeny)

Fast maximum-likelihood phylogenetic tree construction.

**Ubuntu/Debian:**
```bash
sudo apt-get install fasttree
```

**macOS:**
```bash
brew install fasttree
```

**From source:**
```bash
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
sudo mv FastTree /usr/local/bin/
```

### IQ-TREE (Alternative to FastTree)

More accurate phylogenetic inference with bootstrap support.

**Ubuntu/Debian:**
```bash
sudo apt-get install iqtree
```

**macOS:**
```bash
brew install iqtree
```

**From source:**
```bash
wget https://github.com/iqtree/iqtree2/releases/download/v2.2.0/iqtree-2.2.0-Linux.tar.gz
tar -xzf iqtree-2.2.0-Linux.tar.gz
sudo cp iqtree-2.2.0-Linux/bin/iqtree2 /usr/local/bin/
```

### DSSP (Optional)

Secondary structure assignment from PDB files.

**Ubuntu/Debian:**
```bash
sudo apt-get install dssp
```

**macOS:**
```bash
brew install brewsci/bio/dssp
```

### CD-HIT (Optional)

Sequence clustering for redundancy removal.

**Ubuntu/Debian:**
```bash
sudo apt-get install cd-hit
```

**macOS:**
```bash
brew install cd-hit
```

## Verification

Verify your installation:

```python
# Test Python package
import evomotif
print(evomotif.__version__)

# Test all modules
from evomotif import (
    retrieval,
    alignment,
    conservation,
    motif_discovery,
    phylogeny,
    structure,
    variants,
    stats
)
```

Check external tools:

```bash
# Check MAFFT
mafft --version

# Check HMMER
hmmbuild -h
hmmsearch -h

# Check FastTree
FastTree 2>&1 | head -1

# Check IQ-TREE (if installed)
iqtree --version

# Check DSSP (if installed)
mkdssp --version
```

## Troubleshooting

### Python Dependencies

If you encounter issues with Python dependencies:

```bash
# Update pip
pip install --upgrade pip

# Install specific versions
pip install numpy==1.24.0 biopython==1.81

# Use conda for difficult packages
conda install -c conda-forge scipy scikit-learn
```

### External Tool Path Issues

If tools are not found in PATH:

```python
# Specify paths explicitly
from evomotif.alignment import SequenceAligner

aligner = SequenceAligner(mafft_path="/custom/path/to/mafft")
```

### Memory Issues

For large alignments (>2000 sequences):

1. Increase system swap space
2. Use MAFFT auto mode instead of L-INS-i
3. Subsample sequences using CD-HIT

### Permission Issues on Linux

```bash
# Add user to necessary groups
sudo usermod -aG bio $USER

# Make tools executable
chmod +x /path/to/tool
```

## Next Steps

- Read the [Tutorial](tutorials/basic_usage.md)
- Check [Examples](examples/index.md)
- Explore [API Documentation](api.md)
