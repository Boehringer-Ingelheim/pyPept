"""
Tests for the converter module

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


########################################################################################
# Modules
########################################################################################

import unittest

from pyPept.converter import Converter

class TestConverter(unittest.TestCase):
    """
    Class to test some Converter functionalities.

    This should be expanded as time permits.
    """

    def test_constructor(self):
        """
        Method for testing inputs in constructor.
        """
        try:
            Converter(biln="A-C", helm='PEPTIDE1{A.C}$$$$V2.0')
        except Exception as e:
            self.assertIsInstance(e, ValueError)

if __name__ == "__main__":
    unittest.main()
