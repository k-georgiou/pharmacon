"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Module :mod:`pharmacon.constants.formats`.
"""
from typing import Final




__all__ = [
    "supported_trajectory_analysis_topology_formats",
    "supported_structure_analysis_topology_formats",
    "supported_trajectory_analysis_trajectory_formats",
    "supported_pharmacon_ini_formats",
]


supported_trajectory_analysis_topology_formats: Final[frozenset[str]] = frozenset({".parm7", ".prmtop", ".psf", ".dms", ".tpr"})
supported_structure_analysis_topology_formats: Final[frozenset[str]] = frozenset({".mol2", ".sdf", ".smi", ".pdb", ".gro", ".crd"})
supported_trajectory_analysis_trajectory_formats: Final[frozenset[str]] = frozenset({".nc", ".dcd", ".xtc", ".trr"})
supported_pharmacon_ini_formats: Final[frozenset[str]] = frozenset({".ini", ".inp", ".conf"})

