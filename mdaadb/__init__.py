"""
MDAnalysis Analysis Database
A short description of the project.
"""

# Add imports here
from importlib.metadata import version

from .query import Query
from .database import Database, Table, Schema
from .analysis import DBAnalysisManager


__version__ = version("mdaadb-kit")
