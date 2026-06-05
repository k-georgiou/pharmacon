"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

This module provides functionality for working with specific file types
related to Pharmacon, including HDF5 and PTA files.

The module imports `PharmaconHDF5File` and `PharmaconPTAFile` as `PTAFile`.
It also defines the public API through the `__all__` directive.
"""


from .base import PharmaconHDF5File
from .pta import PharmaconPTAFile as PTAFile
from .psa import PharmaconPSAFile as PSAFile


__all__ = [
    "PharmaconHDF5File",
    "PTAFile",
    "PSAFile",
]
