"""
Tests for the conformer module

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
import os

from pyPept.sequence import Sequence
from pyPept.sequence import correct_pdb_atoms
from pyPept.molecule import Molecule
from pyPept.conformer import Conformer
from pyPept.conformer import SecStructPredictor

##########################################################################
# Functions and classes
##########################################################################

def write_molecule(biln):
    """
    Test generation of PDB file
    """

    seq = Sequence(biln)
    seq = correct_pdb_atoms(seq)
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt='ROMol')
    fasta = Conformer.get_peptide(biln)
    ss_input = SecStructPredictor.predict_active_ss(fasta)
    Conformer.generate_conformer(romol, ss_input, generate_pdb=True)
    out_file='structure.pdb'

    with open(out_file, 'r', encoding='utf-8') as file_info:
        data = file_info.read()

    os.unlink(out_file)

    return data

########################################################################################

class TestConformer(unittest.TestCase):
    """
    Class to test some conformer functionalities
    """

    def test_conformer(self):
        """
        Function to check the generation of the PDB files
        """
        self.assertEqual(write_molecule('A-am').split('\n')[0],
                         'ATOM      1  CB  ALA A   1       0.289   0.526  -0.057  \
1.00  0.00           C  ')
        self.assertEqual(write_molecule('A-A-A').split('\n')[0],
                         'ATOM      1  CB  ALA A   1       2.347  -3.516   1.575  \
1.00  0.00           C  ')
        self.assertEqual(write_molecule('C(1,3)-A-A-A-C(1,3)').split('\n')[0],
                         'ATOM      1  N   CYS A   1       0.688  -3.202   4.865  \
1.00  0.00           N  ')
        self.assertEqual(write_molecule('P-E-P-T-I-D-E').split('\n')[0],
                         'ATOM      1  N   PRO A   1       3.980  -0.707   2.723  \
1.00  0.00           N  ')
        self.assertEqual(write_molecule('ac-D-T-H-F-E-I-A-am').split('\n')[0],
                         'ATOM      1  CH3 ACE A   1       6.572  -4.655   2.934  \
1.00  0.00           C  ')
        self.assertEqual(write_molecule('A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)\
-A-A-D(2,3)').split('\n')[0],
                         'ATOM      1  CB  ALA A   1       3.225   1.597  -6.954  \
1.00  0.00           C  ')
        self.assertEqual(write_molecule('A-G-Q-A-A-K(1,3)-E-F-I-A-A.G-L-E-E(1,3)').split('\n')[0],
                         'ATOM      1  CB  ALA A   1      -2.198 -11.275   0.578  \
1.00  0.00           C  ')
        self.assertEqual(write_molecule('N-Iva-F-D-I-meT-N-A-L-W-Y-Aib-K').split('\n')[0],
                         'ATOM      1  ND2 ASN A   1     -12.091   1.882  -4.057  \
1.00  0.00           N  ')


if __name__ == "__main__":
    unittest.main()
