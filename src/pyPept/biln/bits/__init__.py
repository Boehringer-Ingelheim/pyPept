"""
Export of abstract data used in BILN-related code.

From publication: pyPept: a python library to generate atomistic
    2D and 3D representations of peptides
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

from .defaults import BILNConstants
from .abstractions import NumberedMonomer
from .abstractions import BILNSequenceError, BILNMultiError

__all__ = [
    "BILNConstants", "NumberedMonomer",
    "BILNSequenceError", "BILNMultiError",
    ]
