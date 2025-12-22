"""
EvoMotif: Evolution-Driven Framework for Protein Motif Discovery

An open-source pipeline for discovering novel protein motifs through
evolutionary conservation analysis, structural mapping, and phylogenetic inference.

Quick Start:
    >>> import evomotif
    >>> results = evomotif.analyze_protein("p53", "your@email.com")
    >>> print(results.summary())
"""

__version__ = "0.1.0"
__author__ = "Taha"
__license__ = "MIT"

# Import high-level API
from evomotif.pipeline import analyze_protein, EvoMotifPipeline, AnalysisResults

# Import modules for advanced users
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

__all__ = [
    # High-level API (recommended for most users)
    'analyze_protein',
    'EvoMotifPipeline',
    'AnalysisResults',
    # Individual modules (for advanced users)
    'retrieval',
    'alignment',
    'conservation',
    'motif_discovery',
    'phylogeny',
    'structure',
    'variants',
    'stats'
]
