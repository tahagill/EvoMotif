"""
EvoMotif: Evolution-Driven Framework for Protein Motif Discovery

An open-source pipeline for discovering novel protein motifs through
evolutionary conservation analysis, structural mapping, and phylogenetic inference.
"""

__version__ = "0.1.0"
__author__ = "Taha"
__license__ = "MIT"

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
    'retrieval',
    'alignment',
    'conservation',
    'motif_discovery',
    'phylogeny',
    'structure',
    'variants',
    'stats'
]
