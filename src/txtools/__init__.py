"""txtools: transcriptomic, nucleotide-resolution analysis of RNA-seq.

A Python reimplementation of the txtools R package (Garcia-Campos et al.,
Nucleic Acids Research 2024) focused on speed and easy installation.
"""
from .api import *  # noqa: F401,F403
from .api import __all__  # noqa: F401

__version__ = "0.1.0"
