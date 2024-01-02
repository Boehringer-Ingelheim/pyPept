"""
Tests for the molecule module

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

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule

##########################################################################
# Functions and classes
##########################################################################

def try_get_molecule(biln, fmt='Smiles'):
    """
    Test conversion biln -> Molecule
    """

    if biln is not None and biln != 'String':
        seq = Sequence(biln)
    else:
        seq = biln
    mol = Molecule(seq)

    out = mol.get_molecule(fmt=fmt)

    return out

########################################################################################

class TestMolecule(unittest.TestCase):
    """
    Class to test some Molecule functionalities
    """

    def test_get_molecule(self):
        """
        Function to test the generation of SMILES from the molecule
        :return:
        """
        self.assertEqual(try_get_molecule('A-am'), 'C[C@H](N)C(N)=O')
        self.assertEqual(try_get_molecule('A-A-A'), 'C[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O')
        self.assertEqual(try_get_molecule('C(1,3)-A-A-A-C(1,3)'),
                         'C[C@@H]1NC(=O)[C@H](C)NC(=O)[C@@H](N)\
CSSC[C@@H](C(=O)O)NC(=O)[C@H](C)NC1=O')
        self.assertEqual(try_get_molecule('P-E-P-T-I-D-E'),
                         'CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@@H]1CCCN1C(=O)[C@H](CCC(=O)O)NC(=O)\
[C@@H]1CCCN1)[C@@H](C)O)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)O')
        self.assertEqual(try_get_molecule('ac-D-T-H-F-E-I-A-am'), 'CC[C@H](C)[C@H](NC(=O)[C@H](CCC\
(=O)O)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@@H](NC(=O)\
[C@H](CC(=O)O)NC(C)=O)[C@@H](C)O)C(=O)N[C@@H](C)C(N)=O')
        self.assertEqual(try_get_molecule('A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)-A-A-D(2,3)'),
                         'C[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@H]1CC(=O)NCCCC[C@H](N)C(=O)\
N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@H](C(=O)O)CC(=O)NCCCC[C@@H]\
(C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC1=O')
        self.assertEqual(try_get_molecule('A-G-Q-A-A-K(1,3)-E-F-I-A-A.G-L-E-E(1,3)'),
                         'CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CCC(=O)O)\
NC(=O)[C@H](CCCCNC(=O)CC[C@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H]\
(CC(C)C)NC(=O)CN)C(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H]\
(CCC(N)=O)NC(=O)CNC(=O)[C@H](C)N)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O')
        self.assertEqual(try_get_molecule('N-Iva-F-D-I-meT-N-A-L-W-Y-Aib-K'),
                         'CC[C@H](C)[C@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1ccccc1)NC(=O)\
[C@](C)(CC)NC(=O)[C@@H](N)CC(N)=O)C(=O)N(C)[C@H](C(=O)N[C@@H](CC(N)=O)\
C(=O)N[C@@H](C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)\
C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NC(C)(C)C(=O)N[C@@H]\
(CCCCN)C(=O)O)[C@@H](C)O')

if __name__ == "__main__":
    unittest.main()
