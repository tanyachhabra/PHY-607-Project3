"""
ising2d - A package for simulating the 2D Ising model
===================================================

This package provides tools for simulating and analyzing the 
2D Ising model using the Swendsen-Wang cluster algorithm.
"""

from . import model
from . import analysis
from . import visualization
from . import config

__all__ = ['model', 'analysis', 'visualization', 'config']
