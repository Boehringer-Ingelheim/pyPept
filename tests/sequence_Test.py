"""
Unit tests for the sequence module.

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"
__version__ = "1.0"

########################################################################################
# Modules
########################################################################################

import unittest
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule

##########################################################################
# Functions and classes
##########################################################################

def count_atoms(biln):
    """
    Count the number of atoms after generating the molecule.

    :param biln: BILN string of the peptide
    :return: number of atoms in molecule
    """
    seq = Sequence(biln)
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt='ROMol')
    nat = romol.GetNumAtoms()

    return nat

##########################################################################

def sequence_info(biln):
    """
    Get information from the sequence object

    :param biln: biln string of the peptide
    :return: list of monomers and bonds from the molecule

    """

    seq = Sequence(biln)
    # Make sure bonds come back sorted by first residue,
    # first attachment point to allow comparison
    bonds = sorted(seq.s_bonds, key=lambda x: x[1])
    bonds = sorted(bonds, key=lambda x: x[0])

    long_name=[]
    for mon in seq.s_monomers:
        long_name.append(mon['m_name'])

    return long_name, bonds

##########################################################################

class TestSequence(unittest.TestCase):
    """
    Class to test some of the sequence functionalities
    """

    def test_mol_from_sequence(self):
        """
        Test if correct molecule is generated from BILN
        """
        self.assertEqual(count_atoms('A-am'),6)
        self.assertEqual(count_atoms('A-A-A'),16)
        self.assertEqual(count_atoms('C(1,3)-A-A-A-C(1,3)'),28)
        self.assertEqual(count_atoms('P-E-P-T-I-D-E'),56)
        self.assertEqual(count_atoms('ac-D-T-H-F-E-I-A-am'),62)
        self.assertEqual(count_atoms('A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)-A-A-D(2,3)'),74)
        self.assertEqual(count_atoms('A-G-Q-A-A-K(1,3)-E-F-I-A-A.G-L-E-E(1,3)'),106)
        self.assertEqual(count_atoms('N-Iva-F-D-I-meT-N-A-L-W-Y-Aib-K'),113)


    def test_sequence_perception(self):
        """
        Test elements of the sequence object
        """
        self.assertEqual(sequence_info("A-am"),
                        (['Alanine', 'C-Terminal amine'], [[0, 3, 1, 0]]), )
        self.assertEqual(sequence_info('A-A-A'),
                         (['Alanine', 'Alanine', 'Alanine'], [[0, 3, 1, 2], [1, 3, 2, 2]]), )
        self.assertEqual(sequence_info('C(1,3)-A-A-A-C(1,3)'),
                         (['Cysteine', 'Alanine', 'Alanine', 'Alanine', 'Cysteine'],
                          [[0, 3, 4, 3], [0, 4, 1, 2], [1, 3, 2, 2], [2, 3, 3, 2], [3, 3, 4, 0]]), )
        self.assertEqual(sequence_info('P-E-P-T-I-D-E'),
                         (['Proline', 'Glutamic acid', 'Proline', 'Threonine', 'Isoleucine',
                           'Aspartic acid', 'Glutamic acid'],
                          [[0, 5, 1, 0], [1, 6, 2, 0], [2, 5, 3, 4], [3, 5, 4, 5],
                           [4, 6, 5, 0], [5, 5, 6, 0]]), )
        self.assertEqual(sequence_info('ac-D-T-H-F-E-I-A-am'),
                         (['N-Terminal acetic acid', 'Aspartic acid', 'Threonine', 'Histidine',
                           'Phenylalanine','Glutamic acid', 'Isoleucine',
                           'Alanine', 'C-Terminal amine'],
                          [[0, 1, 1, 0], [1, 5, 2, 4], [2, 5, 3, 0], [3, 8, 4, 0], [4, 9, 5, 0],
                           [5, 6, 6, 5], [6, 6, 7, 2], [7, 3, 8, 0]]), )
        self.assertEqual(sequence_info('A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)-A-A-D(2,3)'),
                         (['Alanine', 'Alanine', 'Aspartic acid', 'Alanine', 'Alanine', 'Lysine',
                           'Alanine','Alanine', 'Lysine', 'Alanine', 'Alanine', 'Aspartic acid'],
                          [[0, 3, 1, 2], [1, 3, 2, 0], [2, 3, 8, 6], [2, 5, 3, 2], [3, 3, 4, 2],
                           [4, 3, 5, 0], [5, 6, 11, 3], [5, 7, 6, 2], [6, 3, 7, 2], [8, 7, 9, 2],
                           [9, 3, 10, 2], [10, 3, 11, 0]]), )
        self.assertEqual(sequence_info('A-G-Q-A-A-K(1,3)-E-F-I-A-A.G-L-E-E(1,3)'),
                         (['Alanine', 'Glycine', 'Glutamine', 'Alanine', 'Alanine', 'Lysine',
                           'Glutamic acid','Phenylalanine', 'Isoleucine', 'Alanine', 'Alanine',
                           'Glycine', 'Leucine','Glutamic acid', 'Glutamic acid'],
                          [[0, 3, 1, 0], [1, 2, 2, 6], [2, 7, 3, 2], [3, 3, 4, 2], [4, 3, 5, 0],
                           [5, 6, 14, 4],[5, 7, 6, 0], [6, 6, 7, 0], [7, 9, 8, 5], [8, 6, 9, 2],
                           [9, 3, 10, 2],[11, 2, 12, 5], [12, 6, 13, 0], [13, 6, 14, 0]]), )
        self.assertEqual(sequence_info('N-Iva-F-D-I-meT-N-A-L-W-Y-Aib-K'),
                         (['Asparagine', 'Isovaline', 'Phenylalanine', 'Aspartic acid',
                           'Isoleucine','N-Methyl-Threonine', 'Asparagine', 'Alanine', 'Leucine',
                           'Tryptophan', 'Tyrosine',
                           'Alpha-aminoisobutyric acid (2-aminoalanine)', 'Lysine'],
                          [[0, 6, 1, 4], [1, 5, 2, 0], [2, 9, 3, 0], [3, 5, 4, 5], [4, 6, 5, 4],
                           [5, 6, 6, 5],[6, 6, 7, 2], [7, 3, 8, 5], [8, 6, 9, 0], [9, 12, 10, 7],
                           [10, 8, 11, 3], [11, 4, 12, 0]]), )

    def test_length_operators(self):
        self.assertEqual(Sequence('A-A-A').length(), len(Sequence('A-A-A')))

if __name__ == "__main__":
    unittest.main()
