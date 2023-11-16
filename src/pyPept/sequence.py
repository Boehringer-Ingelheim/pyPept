"""
A class to parse and manipulate BILN sequences.

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
import sys
import copy
import re
import os
import warnings

import string
from collections import defaultdict
from importlib.resources import files

# Third-party libraries
import numpy as np
from rdkit import Chem
from rdkit.Chem import PandasTools


##########################################################################
# Functions and classes
##########################################################################

class SequenceConstants:
    """
    A class to hold defaults values related to pyPept.Sequence objects.
    """
    def_path = "pyPept.data"
    def_lib_filename = "monomers.sdf"
    monomer_join = "-"
    chain_separator = "."
    csv_separator = ","
    helm_polymer = '|'
    max_rgroups = 4


# End of Sequence class-related constants definition.
############################################################

class Sequence:
    """
    This class holds the information for a given sequence in BILN format,
    parsing monomer and bond information.
    """

    ############################################################
    def __init__(self, input_biln, path=SequenceConstants.def_path,
                 monomer_lib=SequenceConstants.def_lib_filename):
        """
        Instantiate a pyPept.Sequence object with the input BILN sequence.

        :param input_biln: BILN representation of the peptide
        :type input_biln: str
        :param path: an override option to specify a new location of a monomer
                     library.
        :type path: str
        :param monomer_lib: an override option to specific a different monomer
                           library file name.
        :type path: str
        """
        # Variables to store the monomers and bonds
        self.s_inputbiln = input_biln
        self.s_mid = -1
        self.s_bonds = []
        self.s_nbonds = -1
        self.s_monomers = []
        self.s_nmonomers = 0
        self.__is_valid = True

        seq = split_outside(self.s_inputbiln,
                            by_element=SequenceConstants.monomer_join,
                            outside='[]')
        seq = SequenceConstants.monomer_join.join(list(filter(len, seq)))
        # Remove dangling '-' around chain breaks
        self.s_biln = re.sub(r'[-]*\.[-]*',
                             SequenceConstants.chain_separator, seq)

        # Read the monomer dictionary
        default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
        monomer_df_filepath = files(path).joinpath(monomer_lib)

        if monomer_df_filepath.is_file() is False:
            monomer_df_filepath = default_monomer_df_filepath

        self.monomer_df = get_monomer_info(str(monomer_df_filepath))

        try:
            # Parse the BILN sequence
            self.__parse_biln()
        except IOError:
            warnings.warn(
                f"Unspecified error in parsing BILN {self.s_inputbiln}")
            sys.exit(1)

        try:
            # Parse the peptide bonds based on the BILN representation
            self.__parse_biln_bonds()
        except IOError:
            warnings.warn(
                f'Exception in parsing BILN bonds. (BILN: {self.s_biln})')
            sys.exit(2)

        # Fix the R groups
        self.__fix_m_rgroups()

    ########################################################################################
    def __parse_biln(self):
        """Split BILN string into individual monomers, do some basic checks,
        then construct the monomers and store the chemical representations
        of the monomers in the Sequence object.
        """

        biln = self.s_biln
        chains = split_outside(biln,
                               by_element=SequenceConstants.chain_separator,
                               outside='[]')

        self.s_chains = {}
        self.s_chains["s_nChains"] = len(chains)
        self.s_chains["s_cType"] = [None] * len(chains)
        self.s_chains["s_monomerIDs"] = [None] * len(chains)

        m_idx = 0
        # Iterate over the chains
        for num_chain, chain in enumerate(chains):
            residues = split_outside(chain,
                                     by_element=SequenceConstants.monomer_join,
                                     outside='[]')
            monomer_types = set()
            monomer_ids = []
            for res in residues:
                # Extract the residue information
                resname = re.sub(r'\(\d+,\d+\)', '', res)
                resname = re.sub(r'^\[(.*)\]$', '\\1', resname)

                # Check if name exists in MonomerDic

                if resname not in self.monomer_df.index:
                    warnings.warn(
                        f"Monomer {resname} in BILN not found in MonomerDic")
                    warnings.warn(
                        "Need to check BILN or update MonomerDic before proceeding.")
                    sys.exit(3)

                # Add additional information of the monomers
                mm_info = self.monomer_df.loc[resname, :].to_dict()
                mm_value = {'m_name': resname,
                            'm_name_in_biln': res,
                            'm_chainID': num_chain}
                mm_value.update(mm_info)

                # Append the monomer in the sequence object
                self.__append_monomer(mm_value)

                monomer_ids.append(m_idx)
                monomer_types.add(mm_info['m_type'])
                m_idx += 1

            # Save the type of monomer
            if len(monomer_types) == 1:
                if monomer_types == {'aa'}:
                    self.s_chains["s_cType"][num_chain] = 'peptide'
                else:
                    self.s_chains["s_cType"][num_chain] = 'chem'
            else:
                self.s_chains["s_cType"][num_chain] = 'mixed'

            self.s_chains["s_monomerIDs"][num_chain] = monomer_ids

            # end of looping over each residue in a chain
        # end of looping over all chains in the peptide

    # end of internal-use method for parsing BILN.

    ########################################################################################
    def __append_monomer(self, monomer):
        """Place a new monomer at the end of the sequence, i.e.
        append the current monomer to the monomer list.

        :param monomer: monomer dictionary of added monomer
        :type monomer: dictionary
        """
        self.s_mid += 1

        # Fields required for the addition of the monomer
        keys_needed = ['m_name', 'm_abbr', 'm_name_in_biln', 'm_type',
                       'm_subtype', 'm_chainID', 'm_Rgroups',
                       'm_RgroupIdx', 'm_attachmentPointIdx', 'm_romol']

        # check if all necessary keys are available in the monomer definition
        for key in keys_needed:
            try:
                monomer[key]
            except KeyError:
                warnings.warn(f'Key {key} missing in monomer description')
                sys.exit(1)

        # Construct the monomer dictionary
        m_dict = {}
        m_dict.update({key: monomer[key] for key in keys_needed})

        # Update class variables
        self.s_nmonomers += 1
        self.s_monomers.append(m_dict)

    ########################################################################################
    def __parse_biln_bonds(self):
        """
        Function to parse bond information from the BILN string
        """

        # Local variable to store bond information
        bond_info = []
        bond_info_helm = []

        residues = split_outside(self.s_biln,
                                 SequenceConstants.chain_separator +
                                 SequenceConstants.monomer_join, '[]')
        nres = len(residues)
        # Iterate over residues
        for num_res, res in enumerate(residues):

            if num_res < nres - 1:
                if self.__get_monomer_prop('m_chainID', num_res) == \
                        self.__get_monomer_prop('m_chainID', num_res + 1):
                    # We have a bond between the two monomers
                    r1_value = self.__get_monomer_prop(
                        'm_attachmentPointIdx', num_res)
                    r2_value = self.__get_monomer_prop(
                        'm_attachmentPointIdx', num_res + 1)

                    bond = (num_res, r1_value[1], num_res + 1, r2_value[0])
                    self.__add_bond(*bond)

            # Find connections with other chains
            match = re.findall(r"\((\d+),(\d+)\)", res)
            for i, m_value in enumerate(match):
                b_idx, rgrp_idx = int(m_value[0]), int(m_value[1])
                attchment_idx = self.__get_monomer_prop('m_attachmentPointIdx',
                                                        num_res)

                bond_info.append([b_idx, num_res, rgrp_idx,
                                  int(attchment_idx[rgrp_idx - 1])])
                bond_info_helm.append([b_idx, num_res, rgrp_idx])

            # For multiple branching
            match = re.findall(r"\((\d+),(\d+);(\d+),(\d+)\)", res)
            for i, m_value in enumerate(match):
                b1_idx, rgrp1_idx, b2_idx, rgrp2_idx = int(m_value[0]), int(
                    m_value[1]), \
                                                       int(m_value[2]), int(
                    m_value[3])
                attchment_idx = self.__get_monomer_prop('m_attachmentPointIdx',
                                                        num_res)

                bond_info.append([b1_idx, num_res, rgrp1_idx,
                                  int(attchment_idx[rgrp1_idx - 1])])
                bond_info.append([b2_idx, num_res, rgrp2_idx,
                                  int(attchment_idx[rgrp2_idx - 1])])
                bond_info_helm.append([b1_idx, num_res, rgrp1_idx])
                bond_info_helm.append([b2_idx, num_res, rgrp2_idx])

        # Now that we have the bond info we collect additional information
        nbonds = int(len(bond_info))

        if (nbonds % 2) == 1:
            raise ValueError("Must be an even number of extra bonds.")

        if nbonds > 0:
            # Convert to numpy array as this way it is trivial to get to the bond identifiers
            array_val = np.asarray(bond_info)
            bond_identifiers = set(array_val[:, 0])

            bond = [None] * int(nbonds / 2)

            # Find matching bond information based on bond identifier (first number in bond info)
            for idx, bond_val in enumerate(bond_identifiers):
                for i in range(nbonds):
                    if bond_info[i][0] == bond_val:
                        try:
                            bond[idx].append(bond_info[i])
                        except:
                            bond[idx] = [bond_info[i]]

            # Collect the data needed to add to bond information in Sequence
            for bond_val, bondx in enumerate(bond):

                if bondx[0][1] > bondx[1][1]:
                    bondx[0], bondx[1] = bondx[1], bondx[0]

                backbone = [bondx[0][1], bondx[0][3], bondx[1][1], bondx[1][3]]
                self.__add_bond(*backbone)

        # Filter unique bonds
        self.__only_unique_bonds()

    ## end of Sequence.__parse_biln_bonds()

    ########################################################################################
    def __fix_m_rgroups(self):
        """
        Use the bond information to replace R-groups where necessary.
        """
        # collect the m_Rgroup information
        all_rgroups = []
        all_attch_idx = []
        for mon in self.s_monomers:
            m_rgroups = copy.copy(mon['m_Rgroups'])
            m_attch_idx = mon['m_attachmentPointIdx']
            all_rgroups.append(m_rgroups)
            all_attch_idx.append(m_attch_idx)

        for bond in self.s_bonds:
            r1_value, a1_value, r2_value, a2_value = bond
            attach1 = all_attch_idx[r1_value]
            rgroup1 = attach1.index(a1_value)
            attach2 = all_attch_idx[r2_value]
            rgroup2 = attach2.index(a2_value)

            # check if N-term
            all_rgroups[r1_value][rgroup1] = None
            all_rgroups[r2_value][rgroup2] = None

        # Now update the m_Rgroups and romol in the sequence object
        for i, mon in enumerate(self.s_monomers):
            mon['m_Rgroups'] = all_rgroups[i]
            m_rgroup_idx = mon['m_RgroupIdx']
            m_romol = mon['m_romol']
            m_romol = self.__remove_rgroups(all_rgroups[i], m_rgroup_idx,
                                            m_romol)
            mon['m_romol'] = m_romol

    ########################################################################################
    def __add_bond(self, m_id1, atom1, m_id2, atom2):
        """
        Add an entry in the bond list

        :param m_id1: index of monomer 1
        :type m_id1: int
        :param m_id2: index of monomer 2
        :type m_id2: int
        :param atom1: index of bond atom in monomer 1
        :type atom1: int
        :param atom2: index of bond atom in monomer 2
        :type atom2: int
        """
        self.s_nbonds += 1
        self.s_bonds.append([m_id1, atom1, m_id2, atom2])

    ########################################################################################
    def get_monomer(self, idx):
        """
        Return the dictionary of ith monomer by index

        :param idx: index of monomer
        :type idx: int
        """

        mon = None
        try:
            mon = self.s_monomers[idx]
        except IndexError:
            warnings.warn(f'Cannot access monomer index {idx} in sequence')

        return mon

    ########################################################################################
    def __get_monomer_prop(self, property_val, idx):
        """
        Return a property of the i-th monomer

        :param property_val: which property to fetch
        :type property_val: str
        :param idx: index of monomer
        :type idx: int
        :return: desired property
        """

        m_prop = None
        try:
            m_dict = self.get_monomer(idx)
        except ValueError:
            m_prop = None

        try:
            m_prop = m_dict[property_val]
        except ValueError:
            warnings.warn(
                f'Cannot access property {property_val} in monomer {idx}')

        return m_prop

    ############################################################################
    def length(self):
        """
        Return the length of the sequence
        """
        return len(self.s_monomers)

    ############################################################################
    def __len__(self):
        """
        Return the length of the sequence.
        """
        return self.length()

    ############################################################################
    def __only_unique_bonds(self):
        """
        In __parseBILNBonds it might happen that a bond is added twice
        (once from the amide section, once if it is)
        in addition set explicitly by (1,2) sets.
        Here, we filter these bonds out

        :return: filtered bond list with unique bonds only
        """
        bonds = self.s_bonds
        bonds = [list(x) for x in set(tuple(x) for x in bonds)]
        self.s_bonds = bonds

    ########################################################################################
    @staticmethod
    def __remove_rgroups(rgroups, rgroup_idx, romol):
        """
        If an Rgroup in monomer definition is not part of a bond, replace it by the appropriate atom
        :param rgroups: type of atoms in the R-groups
        :type rgroups: str
        :param rgroup_idx: indexes of the atoms in the R-groups
        :type rgroup_idx: int
        :param romol: RDKit molecule object
        :type romol: RDKit molecule
        :return:
        """
        emol = Chem.RWMol(romol)
        idx = []

        for i in range(0, SequenceConstants.max_rgroups):
            if rgroups[i] is not None:
                at_idx = rgroup_idx[i]
                at_type = rgroups[i]

                if at_type == 'H':
                    # Simply remove this R-Group
                    idx.append(at_idx)

                elif at_type == 'OH':
                    # generate an oxygen atom
                    oh_group = Chem.MolFromSmiles('O')
                    at_o = oh_group.GetAtomWithIdx(0)
                    # replace Rgroup by O
                    emol.ReplaceAtom(at_idx, at_o)
                else:
                    warnings.warn(f'R-group {at_type} not known ' + \
                                  '- need to update code in Sequence.__removeRgroups.')
                    raise ValueError

        # Remove the atoms to be removed - reverse list order
        for index in sorted(idx, reverse=True):
            emol.RemoveAtom(index)

        return Chem.Mol(emol)

    ############################################################################
    def is_valid(self):
        """Flag if the initialized pyPept.Sequence is valid.
        This includes check of R-groups present in the monomer dictionary,
        and R-group connectivity when explicitly given.

        :return: bool
        """
        return self.__is_valid

    ############################################################################
    def __bool__(self):
        """
        Returns if a pyPept.Sequence object has been constructed with a valid
        BILN.

        :seealso: is_valid

        >>> s = Sequence('P-E-P')
        >>> bool(s)
        True
        >>> s = Sequence('P-nonsense-P')
        >>> bool(s)
        False
        """
        return self.is_valid()

    ## End of the Sequence class declaration.


##########################################################################
# Additional functions
##########################################################################

def get_monomer_info(path):
    """
    Convert a monomer SDF file to a Pandas dataframe object.

    :param path: os path of the monomers.sdf file
    :type path: os path

    :return: monomer dictionary as a dataframe
    """
    sdf_file = path
    df_group = PandasTools.LoadSDF(sdf_file)

    groups = ['m_Rgroups', 'm_RgroupIdx', 'm_attachmentPointIdx']
    for idx in df_group.index:
        for group in groups:
            change = df_group[group][idx].split(SequenceConstants.csv_separator)
            if group == 'm_Rgroups':
                updated_change = [None if v == 'None' else v for v in change]
            else:
                updated_change = [None if v == 'None' else int(v) for v in
                                  change]
            df_group[group][idx] = updated_change
    df_group = df_group.set_index('symbol')
    df_group = df_group.rename(columns={"ROMol": "m_romol"})

    return df_group


########################################################################
def greekify(mol, name_aa):
    """
    Converts the atom names into relative names using the Greek alphabet.

    :param mol: RDKit molecule object
    :type mol: RDKit molecule
    :param name_aa: name of the AA that will be modified
    :type name_aa: str
    """

    # Some local fixes to name correctly the natural amino acids if required
    fix_aa = {'W': {'CE1': 'CE2', 'NE2': 'NE1', 'CZ1': 'CZ3', 'CH': 'CH2',
                    'CD1': 'CD2', 'CD2': 'CD1'},
              'N': {'OD2': 'OD1', 'ND1': 'ND2'},
              'I': {'CD': 'CD1', 'CG1': 'CG2', 'CG2': 'CG1'},
              'T': {'OG2': 'OG1', 'CG1': 'CG2'},
              'P': {'CG2': 'CG', 'CG1': 'CD'}}

    greek = list('ABGDEZHTIKLMNXOPRS')
    greekdex = defaultdict(list)
    ca_atom = get_atom_by_name(mol, 'CA')

    # Recognize if the atom is not part of the backbone
    for atom in mol.GetAtoms():
        is_backbone = (atom.GetPDBResidueInfo() is not None and
                       atom.GetPDBResidueInfo().GetName().strip() in (
                           'LOWER', 'UPPER', 'N', 'CA', 'C', 'H', 'HA', 'O',
                           'OXT'))
        if atom.GetSymbol() != 'H' and atom.GetSymbol()[
            0] != 'R' and not is_backbone:
            n_atom = len(
                Chem.GetShortestPath(mol, ca_atom.GetIdx(), atom.GetIdx())) - 1
            greekdex[n_atom].append(atom)

    # Iterate over the list of atoms and the greek letters
    for k in greekdex:
        if len(greek) <= k:
            pass
        elif len(greekdex[k]) == 0:
            pass
        elif len(greekdex[k]) == 1:
            # Special case
            new_name = f'{greekdex[k][0].GetSymbol()}{greek[k]}'
            if name_aa in fix_aa:
                if new_name in fix_aa[name_aa]:
                    digit = fix_aa[name_aa][new_name][-1]
                    name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]}{digit}'
                else:
                    name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '
            else:
                name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '

            greekdex[k][0].GetPDBResidueInfo().SetName(name)

        elif len(greekdex[k]) < 36:
            list_atom = list(string.digits + string.ascii_uppercase)[1:]
            for i, atom in enumerate(greekdex[k]):

                new_name = f'{atom.GetSymbol()}{greek[k]}{list_atom[i]}'
                if name_aa in fix_aa:
                    if new_name in fix_aa[name_aa]:
                        if name_aa != 'P':
                            digit = fix_aa[name_aa][new_name][-1]
                            name = f'{atom.GetSymbol(): >2}{greek[k]}{digit}'
                        else:
                            digit = fix_aa[name_aa][new_name][-1]
                            name = f'{atom.GetSymbol(): >2}{digit} '
                    else:
                        name = f'{atom.GetSymbol(): >2}{greek[k]}{list_atom[i]}'
                else:
                    name = f'{atom.GetSymbol(): >2}{greek[k]}{list_atom[i]}'

                greekdex[k][i].GetPDBResidueInfo().SetName(name)

        else:
            pass


############################################################

def split_outside(string, by_element, outside, keep_marker=True):
    """
    Splits a string by delimiter only if outside of a given delimiter

    :param string: string to be split
    :param by_element: delimiter(s) by which to be split
    :param outside: only split if outside of this
    :param keep_marker: if True keep the chunk marker, remove otherwise

    :return splitChains: split string as list
    """
    # by can be more than 1 character
    by_element = list(by_element)

    # if outside is only one character (e.g. ', "), double it for start and end
    if len(outside) == 1:
        outside = outside + outside

    # Special character
    grpsep = chr(29)

    out = ''
    inside = False
    for i in string:
        if i == outside[0]:
            if inside:
                if keep_marker:
                    j = i
                else:
                    j = ''
                inside = False
            else:
                inside = True
                if keep_marker:
                    j = i
                else:
                    j = ''
        elif i == outside[1]:
            inside = False
            if keep_marker:
                j = i
            else:
                j = ''
        else:
            if not inside and i in by_element:
                j = grpsep
            else:
                j = i
        out = out + j

    # Do the final split
    split_chains = out.split(grpsep)
    return split_chains


############################################################

def get_atom_by_name(mol, name):
    """
    Get an atom RDKit object based on a particular name

    :param mol: RDKit object with the molecule
    :type mol: RDKit molecule
    :param name: name of the atom
    :type name: str

    :return: the atom RDKit object
    """

    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetName().strip() == name.strip():
            selected_atom = atom

    return selected_atom


########################################################################################

def get_monomer_codes(df_name):
    """
    Get the PDB codes from the monomer dictionary
    :param df_name: monomer dictionary
    :type df_name: dataframe

    :return: dictionary with the monomer PDB codes
    """
    # Get the values from the monomer dictionary
    monomers = {}
    for idx in df_name.index:
        symbol = df_name.at[idx, 'm_abbr']
        pdb_name = df_name.at[idx, 'pdbName']
        monomers[symbol] = pdb_name

    return monomers


############################################################

def correct_pdb_atoms(seq, path=SequenceConstants.def_path,
                      monomer_lib=SequenceConstants.def_lib_filename):
    """
    Pipeline to assign the correct atom names to the pyPept object

    :param seq: pyPept Sequence object
    :type seq: pyPept.sequence
    :return: the modified pyPept Sequence with correct atom names
    """

    # Special case two main N- and C- terminal caps
    names_cap = {'ac': ['CH3', 'C', 'O'], 'am': ['N']}

    # Read the monomer dataframe
    default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
    monomer_df_filepath = files(path).joinpath(monomer_lib)

    if monomer_df_filepath.is_file() is False:
        monomer_df_filepath = default_monomer_df_filepath

    new_df = get_monomer_info(str(monomer_df_filepath))


    # Get monomer codes
    monomers = get_monomer_codes(new_df)

    # Iterate over the monomers
    mm_list = seq.s_monomers
    for i, monomer in enumerate(mm_list):
        mol = monomer['m_romol']
        name = monomer['m_abbr']
        backbone_smiles = Chem.MolFromSmiles('NCC(=O)')
        if i == len(mm_list) - 1:
            backbone_smiles = Chem.MolFromSmiles('NCC(=O)O')
        backbone_atoms = mol.GetSubstructMatches(backbone_smiles)
        aa_flag = 0
        if backbone_atoms:
            bb_index = list(backbone_atoms[0])
            type_mon = new_df.loc[new_df['m_abbr'] == name, 'm_type'].item()
            if type_mon == 'aa':
                aa_flag = 1

        # Iterate over the atoms
        counter = 0
        counter_non = 0
        for j, atom in enumerate(mol.GetAtoms()):
            if aa_flag == 1:
                pos = -1
                if atom.GetIdx() in bb_index:
                    pos = bb_index.index(atom.GetIdx())
                if pos == 0:
                    atomname = ' N  '
                elif pos == 1:
                    atomname = ' CA '
                elif pos == 2:
                    atomname = ' C  '
                elif pos == 3:
                    atomname = ' O  '
                elif pos == 4:
                    atomname = ' OXT'
                else:
                    counter += 1
                    atomname = f' {atom.GetSymbol()}{counter} '
            else:
                # Special case for main capping groups
                if name in ('ac', 'am'):
                    if atom.GetSymbol()[0] != 'R':
                        if names_cap[name][j] == 'CH3':
                            atomname = f' {names_cap[name][j]}'
                        else:
                            atomname = f' {names_cap[name][j]}  '
                else:
                    counter_non += 1
                    if counter_non < 10:
                        atomname = f' {atom.GetSymbol()}{counter_non} '
                    else:
                        atomname = f' {atom.GetSymbol()}{counter_non}'

            # Assign the atom object to the peptide molecule
            info = atom.GetPDBResidueInfo()
            if info is None:
                atom.SetMonomerInfo(Chem.AtomPDBResidueInfo(atomName=atomname,
                                                            serialNumber=atom.GetIdx(),
                                                            residueName=f'{monomers[name]}',
                                                            residueNumber=i + 1,
                                                            chainId="A"))

        # Rename atoms using the greek nomenclature
        if aa_flag == 1:
            greekify(mol, name)

    return seq

    # end of definition of correct_pdb_atoms()

############################################################
# End of sequence.py
############################################################