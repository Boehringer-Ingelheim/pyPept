"""
Export of BILNParser and related objects needed by applications.

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023

Updated 2025
From presentation: Modelling a (New) Modality: Computational Tools for Peptide Design
Certara User Group Meeting, Frankfurt, 2025
"""

################################################################################
# Authorship
################################################################################

__credits__ = ["J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


################################################################################
# Modules
################################################################################

from .implementation import BILNParser

__all__ = [
    "BILNParser",
    ]
