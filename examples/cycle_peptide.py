#!/usr/bin/env python

'''
Example: cycle peptide
Input: BILN sequence
SecStr: Manually provided
Output: PDB structure with 3D coordinates
'''

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"
__version__ = "1.0"

########################################################################################
# Modules
########################################################################################

from src.pyPept.sequence import Sequence
from src.pyPept.sequence import correct_pdb_atoms
from src.pyPept import Molecule
from src.pyPept import Conformer

# RDKit modules
from rdkit import Chem

# Start the Sequence object
biln = "C(1,3)-A-A-A-C(1,3)"
################################################
## With HELM: Call the converter to change from HELM to BILN
#helm = "PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0"
#b = Converter(helm=helm)
#biln = b.get_biln()
################################################

seq = Sequence(biln)
# Correct atom names in the sequence object
seq = correct_pdb_atoms(seq)

# Loop wit the included monomers
mm_list = seq.s_monomers
for i, monomer in enumerate(mm_list):
    mon = monomer['m_romol']

# Generate the RDKit object
mol = Molecule(seq, depiction='local')
romol = mol.get_molecule(fmt='ROMol')
print("The SMILES of the peptide is: {}".format(Chem.MolToSmiles(romol)))

# Create the peptide conformer with corrected atom names and secondary structure
secstruct = "-----"
# Generate the conformer
romol = Conformer.generate_conformer(romol, secstruct, generate_pdb=True)
