"""
Class to create a rdkit molecule from BILN sequence information.

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
# pylint: disable=E1101

# System libraries
import warnings

# Third-party libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdDepictor

##########################################################################
# Functions and classes
##########################################################################
class Molecule:
    """
    Wrapper class around a rdkit ROMol object, with customization for
    peptides.
    """

    ############################################################################
    def __init__(self, sequence=None, depiction='local'):
        """
        Initialize a Molecule object, optionally with a Sequence object.
        :param sequence:  input sequence to be converted to a molecule
        :type sequence: pyPept.Sequence object.
        :param depiction: method to generate a 2D image
                          The local method is used by default.
                          The other value is 'rdkit'
        :type depiction: str , choose from 'rdkit' and 'local'
        """

        self.mol = []
        self.offset = []
        self.bondlist = []
        self.monomers = []
        if depiction not in ('rdkit', 'local'):
            raise ValueError(
                f"Depiction was {depiction}, expected 'rdkit' or 'local'.")
        self.depiction = depiction

        # Main function, will modify data members defined above.
        self.__from_sequence(sequence)

        if not isinstance(self.mol, Chem.rdchem.Mol):
            raise RuntimeError('pyPept.Molecule initialization failure: ' +
                'problem initializing rdkit.ROMol')

    ############################################################################
    def __generate_offset_list(self):
        """
        Generate an offset list - each list entry contains the number of atoms
        in all of the preceding monomers
        """
        mons = self.monomers

        offset = [0]
        for i,val in enumerate(mons):
            mol = val['m_romol']
            n_ats = mol.GetNumAtoms()
            offset.append(n_ats + offset[i])

        self.offset = offset

    ############################################################################
    def __combine_all_monomers(self):
        """
        Combine all monomers in a single molecule object.

        Completes internal initialization to prepare an editable molecule.
        """
        mons = self.monomers

        mol = [0]
        for i,val in enumerate(mons):
            monomer = val['m_romol']
            if i == 0:
                mol = monomer
            else:
                mol = Chem.CombineMols(mol, monomer)

        self.mol = Chem.RWMol(mol)

    ############################################################################
    def __add_bonds_to_mol(self):
        """
        Based on bond information, add bonds between the monomers in the
        sequence.
        """

        offset = self.offset

        for bond in self.bondlist:
            m1_value, at1, m2_value, at2 = bond
            at1 = offset[m1_value] + at1
            at2 = offset[m2_value] + at2

            self.mol.AddBond(at1, at2, Chem.BondType.SINGLE)

    ########################################################################################
    def __fixDihedrals(self):
        """
        Fast fix to get a reasonable 2D representation of the Molecule.

        :return: None
        """

        # Define a number of substructures for phi, psi, amide group and sidechains
        psi_mol = Chem.MolFromSmiles('NCC(=O)N')
        phi_mol = Chem.MolFromSmiles('C(=O)NCC(=O)')
        amid_mol = Chem.MolFromSmiles('CC(=O)NC')
        sc_mol = Chem.MolFromSmarts('[CH][CX4][CH2][#6]')

        psi_matches = self.mol.GetSubstructMatches(psi_mol)
        phi_matches = self.mol.GetSubstructMatches(phi_mol)
        amid_matches = list(self.mol.GetSubstructMatches(amid_mol))
        sidechain_matches = list(self.mol.GetSubstructMatches(sc_mol))
        sidechain_matches1 = list(self.mol.GetSubstructMatches(
            Chem.MolFromSmarts('N[CX4][CX4][#6]')))
        [sidechain_matches.append(sm1) for sm1 in sidechain_matches1]

        # Fix for Proline
        proline = Chem.MolFromSmiles('CC(=O)N1C(C(=O))CCC1')
        proline_matches = list(self.mol.GetSubstructMatches(proline))
        [amid_matches.append(sm1) for sm1 in proline_matches]

        # Make sure all open-chain angles are 120 deg
        any_angle = Chem.MolFromSmarts('[*]~[*x0]~[*]')
        all_angles = list(self.mol.GetSubstructMatches(any_angle))

        # This is for exocyclic bonds (e.g. phenol OH)
        any_angle1 = Chem.MolFromSmarts('[*]~[*]~[*x0]')
        all_angle1 = list(self.mol.GetSubstructMatches(any_angle1))

        [all_angles.append(sm1) for sm1 in all_angle1]

        # Search for atoms that have 4 connection points, then remove these
        #   from the all_angles list
        bond4 = Chem.MolFromSmarts('[*D4]')
        bond4 = list(self.mol.GetSubstructMatches(bond4))
        new_angles = []
        for b4 in bond4:
            b4 = list(b4)[0]

            for idx in range(len(all_angles) - 1, -1, -1):
                angle = all_angles[idx]

                if angle[1] == b4:
                    all_angles.remove(angle)
                    new_angles.append(angle)

        new_angles = set(new_angles)

        # Get conformers
        confs = self.mol.GetConformers()
        if len(confs) != 1:
            warnings.warn(f"{len(confs)} conformers for molecule, expected one.")

        # Update the dihedrals
        count_error = 0
        for conf in confs:
            for pm in psi_matches:
                psi = [pm[0], pm[1], pm[2], pm[4]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *psi, 180.)
                except:
                    count_error+=1

            for pm in phi_matches:
                phi = [pm[0], pm[2], pm[3], pm[4]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *phi, 180.)
                except:
                    count_error+=1

            for pm in amid_matches:
                phi = [pm[0], pm[1], pm[3], pm[4]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *phi, 180.)
                except:
                    count_error+=1

            for pm in sidechain_matches:
                phi = [pm[0], pm[1], pm[2], pm[3]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *phi, 180.)
                except:
                    count_error+=1

            for am in all_angles:
                try:
                    Chem.rdMolTransforms.SetAngleDeg(conf, *am, 120.)
                except:
                    count_error+=1
    # end of Molecule.__fixDihedrals()

    ########################################################################################
    def __from_sequence(self, sequence):
        """
        Function to convert a pyPept.Sequence object into a rdkit mol object.

        :param sequence: Input sequence
        :type sequence: pyPept.Sequence or None

        :return: rdkit.Chem.ROMol object.
        """

        if sequence is None:
            return None
        if not sequence.is_valid():
            warnings.warn('Invalid sequence object given! Returning None.')
            return None

        self.monomers = sequence.s_monomers
        self.bondlist = sequence.s_bonds

        # Step 0: generate a atom offset list: offset[i] = Sum (number of atoms up to monomer i-1)
        self.__generate_offset_list()

        # Step 1: put all monomers into a single molecule
        self.__combine_all_monomers()

        # Step 2: from the bond list, add the bonds into molecule
        self.__add_bonds_to_mol()

        # Step3: remove the Rgroups
        self.mol = Chem.DeleteSubstructs(self.mol, Chem.MolFromSmarts('[#0]'))

        # Step 4: sanitize and generate 2D coords
        Chem.SanitizeMol(self.mol)

        # Compute 2D coordinates
        if self.depiction == 'rdkit':
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(self.mol, clearConfs=True)

        if self.depiction=='local':
            AllChem.Compute2DCoords(self.mol, clearConfs=True)
            Chem.AssignAtomChiralTagsFromStructure(self.mol)
            self.__fixDihedrals()

        return self.mol

    ########################################################################################
    def write_molecule(self, fmt=None, out_file=None):
        """
        Write a molecule to a file or stream in a supported format.

        :param fmt: molecule output format. Allowed formats: SDF, PDB, Smiles, ROMol
                          - SDF: string in the format of a sd-file which can be written to file
                          - PDB: string in the format of a pdb-file
                          - Smiles: smiles representation of sequence
                          - ROMol: rdkit.Chem.ROMol object
        :type fmt: string
        :param out_file: output file to be written
                            if outfile = '-': result is written to stdout.
                            if outfile = None: result is returned to calling context.
        :type out_file: string
        """
        mol = self.get_molecule(fmt=fmt)

        if out_file is None:
            return mol
        if out_file == '-':
            print(mol)
        else:
            with open(out_file, 'w', encoding='utf-8') as out_f:
                out_f.write(mol)

    ########################################################################################
    def get_molecule(self, fmt=None):
        """
        Get a molecule representation of the initialized Sequence object.

        :param fmt: molecule output format. Allowed formats: SDF, PDB, Smiles, ROMol
                          - SDF: string in the format of a sd-file which can be written to file
                          - PDB: string in the format of a pdb-file
                          - Smiles: smiles representation of sequence
                          - ROMol: rdkit.Chem.ROMol object
        :type fmt: string
        """

        if fmt == 'SDF':
            mol = Chem.MolToMolBlock(self.mol, includeStereo=True)
        elif fmt == 'PDB':
            mol = Chem.MolToPDBBlock(self.mol)
        elif fmt == 'Smiles':
            mol = Chem.MolToSmiles(self.mol)
        elif fmt == 'ROMol':
            mol = self.mol
        else:
            warnings.warn('fmt must be one of ' +
                '"SDF", "PDB", "Smiles", "ROMol" but was ' +
                f'{fmt}. Return value will be None.')
            mol = None

        return mol

    # End of the Molecule class declaration.

############################################################
# End of molecule.py
############################################################
