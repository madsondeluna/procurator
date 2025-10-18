"""
Procurator: A bioinformatics pipeline for protein analysis and annotation.
"""

__version__ = "0.1.0"

from .main import main
from .io import read_fasta, write_fasta
from .protein import Protein
from .analysis import analyze_sequence
from .orf_finding import find_orfs
from .annotation import annotate_sequence

__all__ = [
    "main",
    "read_fasta",
    "write_fasta",
    "Protein",
    "analyze_sequence",
    "find_orfs",
    "annotate_sequence",
]
