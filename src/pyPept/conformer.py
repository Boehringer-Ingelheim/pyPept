"""
Module to predict peptide conformer using pyPept functions

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""

#############################################################################
# Authorship
#############################################################################

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


#############################################################################
# Modules
#############################################################################
# pylint: disable=E1101

# System and Third-party libraries
import re
import warnings
import math
import os
import sys
from importlib.resources import files

# RDKit modules
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DistanceGeometry
from rdkit.Chem import rdDistGeom

# Biopython
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

# pyPept modules
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from pyPept.sequence import SequenceConstants
from pyPept.sequence import get_monomer_info

##########################################################################
# Functions and classes
##########################################################################

class ConformerConstants:
    """
    A class to hold defaults values related to pyPept.Conformer objects.
    """
    def_path = "pyPept.data"
    def_matrix_filename = "matrix.txt"
    def_ss_filename = 'total_SS.txt'
    aminoacids = {"ALA": "A", "ASP": "D", "GLU": "E", "PHE": "F", "HIS": "H",
                  "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "GLY": "G",
                  "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
                  "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y", "CYS": "C"}
    coil_symbol = "C"
    helix_symbol = "H"
    sheet_symbol = "E"
    blank_symbol = "-"
    expected_ss_symbols = (
        coil_symbol, helix_symbol, sheet_symbol,
        'G', 'T', 'S', 'I', 'B',
        blank_symbol)

# End of Conformer class-related constants definition.
##########################################################################

class Conformer:
    """
    Class with static methods to complement the prediction of conformers
    """

    ########################################################################################
    @staticmethod
    def get_peptide(biln, path=SequenceConstants.def_path,
                          monomer_lib=SequenceConstants.def_lib_filename):
        """
        Obtain a FASTA version of the peptide. 
        Unknown NNAAs are replaced by alanine monomers.

        :param biln: BILN representation of the peptide
        :type biln: str
        :param path: Path of the data folder (if called externally)
        :type path: str (class constant)
        :param monomer_lib: name of the monomer file
        :type monomer_lib: str (class constant)

        :return: A fasta sequence of a natural peptide counterpart - str
        """
        aa_dict = {'G': 'NCC(=O)O',
                  'A': 'N[C@@]([H])(C)C(=O)O',
                  'R': 'N[C@@]([H])(CCCNC(=N)N)C(=O)O',
                  'N': 'N[C@@]([H])(CC(=O)N)C(=O)O',
                  'D': 'N[C@@]([H])(CC(=O)O)C(=O)O',
                  'C': 'N[C@@]([H])(CS)C(=O)O',
                  'E': 'N[C@@]([H])(CCC(=O)O)C(=O)O',
                  'Q': 'N[C@@]([H])(CCC(=O)N)C(=O)O',
                  'H': 'N[C@@]([H])(CC1=CN=C-N1)C(=O)O',
                  'I': 'N[C@@]([H])(C(CC)C)C(=O)O',
                  'L': 'N[C@@]([H])(CC(C)C)C(=O)O',
                  'K': 'N[C@@]([H])(CCCCN)C(=O)O',
                  'M': 'N[C@@]([H])(CCSC)C(=O)O',
                  'F': 'N[C@@]([H])(Cc1ccccc1)C(=O)O',
                  'P': 'N1[C@@]([H])(CCC1)C(=O)O',
                  'S': 'N[C@@]([H])(CO)C(=O)O',
                  'T': 'N[C@@]([H])(C(O)C)C(=O)O',
                  'W': 'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O',
                  'Y': 'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O',
                  'V': 'N[C@@]([H])(C(C)C)C(=O)O'}

        # Read the monomer dataframe
        default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
        monomer_df_filepath = files(path).joinpath(monomer_lib)

        if monomer_df_filepath.is_file() is False:
            monomer_df_filepath = default_monomer_df_filepath

        new_df = get_monomer_info(str(monomer_df_filepath))

        try:
            if SequenceConstants.chain_separator in biln:
                m_seq = biln.split(SequenceConstants.chain_separator)[0]
            else:
                m_seq = biln
        except ValueError:
            warnings.warn(f"No main peptide was detected for peptide \
                          with BILN: {biln}")
            sys.exit(1)

        # Loop through the list of monomers of the main peptide
        total_monomers = []
        monomers = m_seq.split(SequenceConstants.monomer_join)
        for res in monomers:
            mon = re.sub(r'\(\d+,\d+\)', '', res)
            if mon in aa_dict:
                total_monomers.append(mon)
            else:
                # Check if a natural analog is available in the dataframe
                type_mon = new_df.loc[new_df['m_abbr'] == mon, 'm_type'].item()
                if type_mon != 'cap':
                    nat_analog = new_df.loc[new_df['m_abbr'] == mon, 'natAnalog'].item()
                    if nat_analog != "X":
                        total_monomers.append(nat_analog)
                    else:
                        # Run a similarity calculation to chek the most similar AAs
                        biln = f'{mon}'
                        seq = Sequence(biln)
                        mol = Molecule(seq)
                        romol = mol.get_molecule(fmt='ROMol')
                        mol1 = Chem.RemoveHs(romol)

                        sim_total={}
                        for aa_value in aa_dict:
                            mol2 = Chem.MolFromSmiles(aa_dict[aa_value])
                            fp1 = FingerprintMols.FingerprintMol(mol1)
                            fp2 = FingerprintMols.FingerprintMol(mol2)
                            smiles_similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                            sim_total[aa_value] = smiles_similarity

                        temp_order = sorted(sim_total.items(), key=lambda x: x[1])
                        aa_value = temp_order[-1][0]
                        val = float(temp_order[-1][1])
                        if val >= 0.5:
                            total_monomers.append(aa_value)
                        else:
                            # Assign an alanine if no similar AAs are found
                            total_monomers.append("A")
                            warnings.warn(f"Because no natural amino acids were \
                            mapped, an alanine was assigned instead \
                            of {mon}")

        fasta = ''.join(total_monomers)
        return fasta

    ########################################################################################
    @staticmethod
    def fix_hydrogen_atom_names(romol):
        """
        Function to rename hydrogen atoms according to the PDB defaults

        :param romol: RDKit molecule object
        :type romol: RDKit molecule
        """

        # Check the number of hydrogens at each heavy atom
        num_hydrogens = [None] * romol.GetNumAtoms()
        count_hydrogens = [None] * romol.GetNumAtoms()

        # For each heavy atom count the number of hydrogens present and store in numH
        for atom in romol.GetAtoms():

            if atom.GetSymbol() != 'H':
                idx = atom.GetIdx()
                name = atom.GetPDBResidueInfo().GetName()
                i_h = atom.GetTotalNumHs(includeNeighbors=True)
                num_hydrogens[idx] = i_h

                if i_h >= 3:
                    count_hydrogens[idx] = 0
                else:
                    count_hydrogens[idx] = 1
                    if name.strip() in ['NH1', 'NH2']:
                        count_hydrogens[idx] = 0

        # For each hydrogen, get the name of the atom and modify its own name
        for atom in romol.GetAtoms():
            if atom.GetSymbol() == "H":
                bond = atom.GetBonds()[0]
                if bond.GetBeginAtom().GetIdx() == atom.GetIdx():
                    heavy = bond.GetEndAtom()
                else:
                    heavy = bond.GetBeginAtom()

                resinfo = heavy.GetPDBResidueInfo()
                atom.SetMonomerInfo(resinfo)
                heavyname_full = heavy.GetPDBResidueInfo().GetName()
                heavyname = heavyname_full.strip()[1:]

                idx = heavy.GetIdx()
                if num_hydrogens[idx] > 1:
                    count_hydrogens[idx] = count_hydrogens[idx] + 1
                    number = str(count_hydrogens[idx])
                else:
                    if heavyname_full.strip() in ['NH1', 'NH2']:
                        number = "1"
                    else:
                        number = ""

                atomname = 'H' + heavyname + number
                atomname = f'{atomname:>4}'
                atom.GetPDBResidueInfo().SetName(atomname)

    ########################################################################################
    @staticmethod
    def generate_conformer(romol, ss_value, generate_pdb=False,
                           output_name='structure'):
        """
        Generate the conformer with SS restraints and with the correct atom naming

        :param romol: pyPept and RDKit molecular object
        :type romol: pyPept.molecule
        :param ss_value: SS predicted or provided by the user
        :type ss_value: str
        :param generate_pdb: Flag to generate or not a PDB file
        :param generatePDB: bool
        :param output_name: Name of the PDB file generated. A default value is provided
        :type output_name: str

        :return: RDKit mol object with the conformer following SS restraints and correct atom names
        """
        
        # Sanity check on SS symbols:
        invalid = [v for v in ss_value 
                    if v not in ConformerConstants.expected_ss_symbols]
        if len(invalid) > 0:
            raise RuntimeError(
                "%i invalid secondary structure symbols (not %s): %s ." %
                (len(invalid), 
                 " ".join(ConformerConstants.expected_ss_symbols),
                 " ".join(invalid)))

        # Fix hydrogen atom names
        romol = Chem.AddHs(romol)
        Conformer.fix_hydrogen_atom_names(romol)

        # Get the backbone atoms to assign the SS restraints
        backbone_smiles = Chem.MolFromSmiles('NCC(=O)')
        backbone_atoms = romol.GetSubstructMatches(backbone_smiles)
        bounds = rdDistGeom.GetMoleculeBoundsMatrix(romol)

        # Assignment of helical conformation
        for i, ele in enumerate(ss_value):
            if i < len(ss_value) - 4:
                if ele == 'H':  # and ss[i + 4] == 'H':
                    ind_o = backbone_atoms[i][-1]
                    ind_n = backbone_atoms[i + 4][0]
                    bounds[ind_o, ind_n] = 3.2

        # Assignment of beta sheets
        segments = []
        flag = 0
        for i, ele in enumerate(ss_value):
            if i == 0:
                fragment = []
            if ele == 'E':
                flag = 1
                fragment.append(i)
            else:
                if flag == 1:
                    segments.append(fragment)
                    flag = 0
                fragment = []

        if len(segments) == 2 and len(segments[0]) == len(segments[1]):
            for i, ele in enumerate(segments[0]):
                ind_o = backbone_atoms[ele][-1]
                ind_n = backbone_atoms[segments[1][(i + 1) * -1]][0]
                bounds[ind_o, ind_n] = 3.2

        # Generate the new distance matrix and predict the conformer
        try:
            DistanceGeometry.DoTriangleSmoothing(bounds)
            parameters = AllChem.ETKDGv3()
            parameters.randomSeed = 0xf00d
            parameters.SetBoundsMat(bounds)
            parameters.useRandomCoords = True
            AllChem.EmbedMolecule(romol, parameters)
            AllChem.UFFOptimizeMolecule(romol)
            pdb_mol = Chem.MolToPDBBlock(romol)
        except FileExistsError:
            warnings.warn("Failed to generate the conformer in RDKit")
            sys.exit(1)

        # Store the conformer in a new pdb file
        if generate_pdb:
            # Saving the RDKit PDB file to correct the hydrogens later with BioPython
            pdb_predict = open(f'{output_name}.pdb', 'w', encoding="utf8")
            pdb_predict.write(pdb_mol)
            pdb_predict.close()

            # Fix the hydrogen positions in the final PDB file
            parser = PDBParser()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                reference = parser.get_structure('REF', f'{output_name}.pdb')

            chain = reference[0]['A']

            ids = []
            for atom in chain[len(ss_value)]:
                atom_id = atom.id
                if atom_id == "HXT":
                    ids.append(atom_id)
            for i in ids:
                chain[len(ss_value)].detach_child(i)

            # Saving the new structure
            io_output = PDBIO()
            io_output.set_structure(reference)
            io_output.save(f"{output_name}.pdb")

        return romol

    ## End of the Conformer class declaration.

########################################################################################
class SecStructPredictor:
    """
    Class with static methods to predict the secondary structure
    """

    ########################################################################################
    @staticmethod
    def get_matrix(path):
        """
        Read the matrix file to run the similarity analysis
        :param path: os path of the matrix file
        :type path: os path

        :return: matrix dictionary
        """
        matrix_temp = {}
        with open(path, 'r') as file:
            for line in file:
                fields = line.split()
                key = (fields[0], fields[1])
                matrix_temp[key] = int(fields[2])

        return matrix_temp

    ########################################################################################
    @staticmethod
    def align_position_matrix(peptide1, peptide2,
                              path=ConformerConstants.def_path,
                              matrix_lib=ConformerConstants.def_matrix_filename):
        """
        Align the peptides based on a matrix

        :param peptide1
        :type peptide1: str
        :param peptide2
        :type peptide2: str
        :param path: Path of the data folder (if called externally)
        :type path: str (class constant)
        :param matrix_lib: name of the matrix file
        :type matrix_lib: str (class constant)

        :return: score of the alignment
        """


        # Read the matrix file
        default_matrix_lib_filepath = files(ConformerConstants.def_path).joinpath(ConformerConstants.def_matrix_filename)
        matrix_lib_filepath = files(path).joinpath(matrix_lib)

        if matrix_lib_filepath.is_file() is False:
            matrix_lib_filepath = default_matrix_lib_filepath

        matrix = SecStructPredictor.get_matrix(str(matrix_lib_filepath))

        score_matrix = 0

        for i, pep_seq in enumerate(peptide1):
            # Generate the tuples with the pair of amino acids to compare
            pair1 = (pep_seq, peptide2[i])
            pair2 = (peptide2[i], pep_seq)

            # Obtain the score from the matrix and sum up the values
            if pair1 in matrix:
                value = matrix[pair1]
            else:
                value = matrix[pair2]

            score_matrix += value

        return score_matrix

    ########################################################################################
    @staticmethod
    def similarity_pair(peptide1, peptide2, path = ConformerConstants.def_path,
                        matrix_lib = ConformerConstants.def_matrix_filename):
        """
        Function to calculate similarity between two peptide sequences

        :param peptide1
        :type peptide1: str
        :param peptide2
        :type peptide2: str
        :param path: Path of the data folder (if called externally)
        :type path: str (class constant)
        :param matrix_lib: name of the matrix file
        :type matrix_lib: str (class constant)

        :return: Similarity value - float
        """

        # Alignment between peptide1 and peptide2
        score1_2 = SecStructPredictor.align_position_matrix(peptide1, peptide2,
                                            path=path, matrix_lib=matrix_lib)
        score1_1 = SecStructPredictor.align_position_matrix(peptide1, peptide1,
                                            path=path, matrix_lib=matrix_lib)
        score2_2 = SecStructPredictor.align_position_matrix(peptide2, peptide2,
                                            path=path, matrix_lib=matrix_lib)

        # Calculate similarity value
        sim_val = float(score1_2) / math.sqrt(float(score1_1 * score2_2))

        # Return similarity
        return sim_val

    ############################################################################
    @staticmethod
    def predict_active_ss(sequence, threshold=0.6, frag_size=5,
                          path = ConformerConstants.def_path,
                          matrix_lib = ConformerConstants.def_matrix_filename,
                          ss_lib = ConformerConstants.def_ss_filename):
        """
        Predict Secondary Structure by similarity with active peptides from PDB

        :param sequence: Peptide sequence
        :type sequence: str
        :param threshold: Similarity threshold. Default of 0.6
        :type threshold: float
        :param frag_size: Fragment size to do the comparison. Default is 5
        :type frag_size: int
        :param path: Path of the data folder (if called externally)
        :type path: str (class constant)
        :param matrix_lib: name of the matrix file
        :type matrix_lib: str (class constant)
        :param ss_lib: name of the secondary structure file
        :type ss_lib: str (class constant)

        :return: The predicted Secondary Structure in a string format - str
        """

        # Read the file with the active peptides Secondary Structure info
        default_ss_filepath = files(ConformerConstants.def_path).joinpath(ConformerConstants.def_ss_filename)
        ss_filepath = files(path).joinpath(ss_lib)

        if ss_filepath.is_file() is False:
            ss_filepath = default_ss_filepath

        total_ss = [x.strip() for x in open(str(ss_filepath),encoding="utf8")]

        if len(sequence) >= 5:
            # Store and filter based on peptides bigger than a particular size
            peptides = {}
            for line in total_ss:
                info = line.split()
                seq_ref = info[0]
                ss_value = info[1]
                if 'NH2' not in seq_ref and 'ACE' not in seq_ref:
                    aa_list = seq_ref.split(SequenceConstants.chain_separator)
                    peptide = ''
                    for aa_name in aa_list:
                        peptide += ConformerConstants.aminoacids[aa_name]
                    if len(peptide) == len(ss_value):
                        if len(peptide) >= frag_size:
                            peptides[peptide] = ss_value

            # Store the Secondary Structure profiles from similar peptides on the dataset
            profiles = []
            # Loop to match all possible fragments
            for seq in peptides:
                if len(seq) > len(sequence):
                    for i in range(0, len(seq) - len(sequence) + 1):
                        fragment = seq[i:i + len(sequence)]
                        score = 0
                        for count, e_frag in enumerate(sequence):
                            if e_frag == fragment[count]:
                                score += 1
                        if score >= int(len(sequence) * 2 / 5):
                            similarity = SecStructPredictor.similarity_pair(
                                sequence, fragment, path=path,
                                                         matrix_lib=matrix_lib)
                            if similarity >= threshold:
                                profiles.append(peptides[seq][i:i + len(sequence)])
                                break
                elif len(seq) == len(sequence):
                    similarity = SecStructPredictor.similarity_pair(
                        sequence, seq, path=path, matrix_lib=matrix_lib)
                    if similarity >= threshold:
                        profiles.append(peptides[seq])
                elif len(seq) < len(sequence):
                    for i in range(0, len(sequence) - len(seq) + 1):
                        fragment = sequence[i:i + len(seq)]
                        score = 0
                        for count, e_frag in enumerate(seq):
                            if e_frag == fragment[count]:
                                score += 1
                        if score >= int(len(seq) * 2 / 5):
                            similarity = SecStructPredictor.similarity_pair(
                                seq, fragment, path=path, matrix_lib=matrix_lib)
                            if similarity >= threshold:
                                start = ConformerConstants.coil_symbol * i
                                end = ConformerConstants.coil_symbol * (
                                    len(sequence) - len(seq) - i)
                                profiles.append(start + peptides[seq] + end)
                                break

            # Final assignment based on pattern frequency per position
            final_ss = ''
            if profiles:
                for i, aa_name in enumerate(sequence):
                    categories = {}
                    for profile in profiles:
                        cat = profile[i]
                        if cat not in categories:
                            categories[cat] = 1
                        else:
                            categories[cat] += 1
                    marklist = sorted(categories.items(), key=lambda x: x[1])
                    if i in (0, len(sequence) - 1):
                        final_ss += ConformerConstants.blank_symbol
                    else:
                        final_ss += marklist[-1][0]
            else:
                final_ss = ConformerConstants.blank_symbol * len(sequence)
        else:
            final_ss = ConformerConstants.blank_symbol * len(sequence)

        print(f"Predicted Secondary Structure: {final_ss} for main chain: {sequence}")

        return final_ss
    # end of Conformer.predict_active_ss()
# End of SecStructPredictor class.

############################################################
# End of conformer.py
############################################################
