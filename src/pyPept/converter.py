"""
Class to to manipulate BILN and convert to HELM (and vice versa)

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
# pylint: disable=E1101

# System libraries
import copy
import re
import warnings

# pyPept functions
from pyPept.sequence import SequenceConstants
from pyPept.sequence import split_outside

##########################################################################
# Main class
##########################################################################
class Converter:
    """
    Class to convert between BILN and HELM (and viceversa)
    """
    def __init__(self, biln=None, helm=None):
        """Construct a BILN or HELM molecule. They can be read from a file

        :param biln: The BILN format that will be converted to HELM
        :type biln: str
        :param helm: The HELM format that will be converted to BILN
        :type helm: str
        """

        self.polymerinfo = {"chains": [], "bonds": []}

        if biln is not None and helm is not None:
            raise ValueError(
                "Converter cannot be initialized with both BILN and HELM.")
        elif biln is not None:
            self.eval_biln(biln=biln)
        elif helm is not None:
            self.eval_helm(helm=helm)

    ############################################################
    @staticmethod
    def __remove_brackets(sequence):
        """Remove brackets around individual residues

        :param sequence: the peptide sequence
        :type sequence: list

        :return newseq: new list without the brackets
        :type newseq: list
        """
        newseq = []
        for residue in sequence:
            pat = re.sub(r'\[(.*)\]', r'\1', residue)  # remove brackets if necessary
            newseq.append(pat)

        return newseq

    ############################################################
    @staticmethod
    def __split_helm(helm):
        """split HELM string into parts

        :param helm: the helm peptide sequence
        :type sequence: string

        :return helm_parts: a list containing relevant parts from helm
        :type helm_parts: list
        """

        helm_parts = []
        while len(helm):
            match = re.search(r'\$[^,;]', helm)  # this excludes any $ signs in the cxsmiles
            if match is None:
                match = re.search(r'\$', helm)  # search for all other $ signs

            if match is None:
                helm_parts.append(helm)  # done, no more parts
                break

            helm_parts.append(helm[:match.span()[0]])
            helm = helm[match.span()[0] + 1:]  # remove the current part from helm

        if len(helm_parts) == 4:
            helm_parts.append('')

        # helm part 1 == list of simple polymers - split further
        list_of_simple_polymers = helm_parts[0]
        if SequenceConstants.helm_polymer in list_of_simple_polymers:
            list_of_simple_polymers = list_of_simple_polymers.split(
                SequenceConstants.helm_polymer)
        else:
            list_of_simple_polymers = [list_of_simple_polymers]
        helm_parts[0] = list_of_simple_polymers

        # split the bond information into individual chunks
        list_of_connections = helm_parts[1]
        if list_of_connections != '':
            if SequenceConstants.helm_polymer in list_of_connections:
                list_of_connections = list_of_connections.split(
                    SequenceConstants.helm_polymer)
            else:
                list_of_connections = [list_of_connections]
        else:
            list_of_connections = []

        helm_parts[1] = list_of_connections

        return helm_parts

    ############################################################
    def __to_biln(self):
        """Generate BILN from polymerInfo

        :return biln: the generated biln from polymer dictionary
        :type biln: string
        """

        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        # move bond info into monomers
        for ibond, bond in enumerate(bonds):
            c1_value, res1, rgroup1, c2_value, res2, rgroup2 = bond
            chains[c1_value][res1] = f"{chains[c1_value][res1]}({ibond + 1},{rgroup1})"
            chains[c2_value][res2] = f"{chains[c2_value][res2]}({ibond + 1},{rgroup2})"

        # chain everything together
        list_ofsimple_polymers = []
        for chain in chains:
            poly = SequenceConstants.monomer_join.join(chain)
            list_ofsimple_polymers.append(poly)

        biln = SequenceConstants.chain_separator.join(list_ofsimple_polymers)

        if len(biln) == 0:
            biln = None

        return biln

    ############################################################
    def __to_helm(self):
        """Generate a HELM from an internal PolymerInfo dictionary.

        :return helm: the generated helm v2.0 from polymer dictionary
        :type helm: string
        """

        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        # brackets around residues if necessary
        for count_chain, chain in enumerate(chains):
            for count_res, res in enumerate(chain):
                if len(res) > 1:
                    chains[count_chain][count_res] = f"[{res}]"

        # generate HELM elements
        # compile chains
        list_of_simple_polymers = []
        for count_chain, chain in enumerate(chains):
            poly = ".".join(chain)
            list_of_simple_polymers.append(f'PEPTIDE{count_chain + 1}{{{poly}}}')

        # compile bondinfo
        list_of_connections = []
        for bond in bonds:
            c1_val, r1_val, g1_val, c2_val, r2_val, g2_val = bond
            bond_info = f'PEPTIDE{c1_val + 1},PEPTIDE{c2_val + 1},{r1_val + 1}:R{g1_val}-{r2_val + 1}:R{g2_val}'
            list_of_connections.append(bond_info)

        list_of_connections = SequenceConstants.helm_polymer.join(list_of_connections)
        list_of_simple_polymers = SequenceConstants.helm_polymer.join(list_of_simple_polymers)

        if not list_of_simple_polymers:
            helm = None
        else:
            helm = f"{list_of_simple_polymers}${list_of_connections}$$$V2.0"

        return helm

    ############################################################
    def eval_helm(self, helm):
        """Generate internal PolymerInfo dictionary from a HELM string.

        This method fills the polymerinfo dictionary:
        polymerinfo = {"chains": [chain1, chain2, ...],
                       "bonds": [b1, b2, b3]}
        where chains contains the monomer abbreviation for each chain
        and bonds contains the "extra" bonds explicitly given in HELM or BILN as a list

        :param helm: the peptide helm string
        :type helm: string
        """

        # split the helm string into its individual parts
        try:
            list_of_simple_polymers, list_of_connections, groups, annotations,\
            version = self.__split_helm(helm)

        except (ValueError, IndexError):
            warnings.warn(f'problem with HELM string - not enough sections: {helm}')
            warnings.warn(f'need 5, have {len(self.__split_helm(helm))}')
            return None

        if len(list_of_simple_polymers) == 0:
            warnings.warn(f'No simple polymers in HELM string {helm}')
            return None

        # go through the list of simple polymers (first component of HELM string),
        # parse them and put them into polymerinfo["chains"]

        pattern = re.compile(r'{.*}')

        id_val = []
        polymer = []

        for idx, chain in enumerate(list_of_simple_polymers):

            # remove any whitespace
            chain = chain.strip()
            # split each polymer into name/identifier and sequence
            match = pattern.search(chain)
            if match is None:
                warnings.warn('No sequence information found in simple polymer - check HELM')
                warnings.warn(f'Input: {chain}')
                return None

            # split each polymer into name/identifier and sequence
            seq = match.span()
            id_chain = chain[:seq[0]]  # identifier
            if id_chain[0:4] == 'CHEM':
                warnings.warn('polymer contains an explicit chemical entity - \
                missing or not recognized monomer:')
                warnings.warn(chain)
                return None
            if id_chain[0:7] != 'PEPTIDE':
                warnings.warn('non-peptide chains in HELM - probably missing monomer')
                warnings.warn(f'found: {id_chain}')
                return None
            id_chain = int(re.sub('PEPTIDE', '', id_chain))

            poly = chain[seq[0] + 1:seq[1] - 1]  # sequence

            if not poly:
                warnings.warn(f'simple polymer {poly} is declared but not defined - has no length')
                return None

            # split into individual monomers.
            poly = split_outside(poly, SequenceConstants.chain_separator, '[]')

            # remove brackets around monomer abbreviations with more than 1 letter
            poly = self.__remove_brackets(poly)

            id_val.append(id_chain)
            polymer.append(poly)

        self.polymerinfo["chains"] = polymer

        # parse the bond information

        # bond information is of the type: 'PEPTIDE1,PEPTIDE2,1:R1-4:R3'
        # split into individual parts 'id1, id2, res1:rgroup1-res2:rgroup2'
        bonds = []
        if len(list_of_connections) > 0:
            for idx, conn in enumerate(list_of_connections):
                id1, id2, bond = conn.split(',')

                res1, rgroup1, res2, rgroup2 = re.split(r'[-:]', bond)

                # need some reformatting
                id1 = int(id1.replace('PEPTIDE', ''))
                id2 = int(id2.replace('PEPTIDE', ''))

                # translate back to current order of chains in polymerinfo["chains"]
                id1 = id_val.index(id1)
                id2 = id_val.index(id2)

                # start counting of residues at 0
                res1 = int(res1) - 1
                res2 = int(res2) - 1

                # get the number of the Rgroup (keep numbering 1-3)
                rgroup1 = int(rgroup1.replace('R', ''))
                rgroup2 = int(rgroup2.replace('R', ''))

                bonds.append([id1, res1, rgroup1, id2, res2, rgroup2])

        self.polymerinfo["bonds"] = bonds

    ############################################################
    def eval_biln(self, biln):
        """Generate internal PolymerInfo dictionary from a BILN string.

        This method fills the PolymerInfo dictionary:
        polymerinfo = {"chains": [chain1, chain2, ...],
                       "bonds": [b1, b2, b3]}
        where chains contains the monomer abbreviation for each chain
        and bonds contains the "extra" bonds explicitly given in HELM or 
        BILN as a list.

        :param biln: the peptide BILN string
        :type biln: string

        :return: None
        """

        # split BILN into chains
        chains = biln.split(".")

        bondinfo = {}
        list_of_simple_polymers = []

        for cidx, chain in enumerate(chains):
            residues = chain.split(SequenceConstants.monomer_join)

            # go through residues and collect the bond information needed
            polymer = []
            for ridx, res in enumerate(residues):
                match = re.findall(r"\((\d+),(\d+)\)", res)
                resname = re.sub(r"\(.*\)", "", res)
                if match:
                    for m_val in match:
                        bidx = m_val[0]
                        gidx = int(m_val[1])
                        try:
                            bondinfo[bidx].append(cidx)
                        except:
                            bondinfo[bidx] = [cidx]

                        bondinfo[bidx].append(ridx)
                        bondinfo[bidx].append(gidx)

                polymer.append(resname)
            list_of_simple_polymers.append(polymer)

        self.polymerinfo["chains"] = list_of_simple_polymers

        # compile bondinfo
        list_of_connections = []
        for key in bondinfo:
            if len(bondinfo[key]) != 6:
                warnings.warn(f'Error in bond information of BILN: bond {key} not correct')
                return None

            c1_val, r1_val, g1_val, c2_val, r2_val, g2_val = bondinfo[key]

            bond = [c1_val, r1_val, g1_val, c2_val, r2_val, g2_val]
            list_of_connections.append(bond)

        self.polymerinfo["bonds"] = list_of_connections

    ############################################################
    def get_helm(self):
        """Return a HELM string from the PolymerInfo generated at
        instantiation.
        
        :return: str"""
        helm = self.__to_helm()
        return helm

    ############################################################
    def get_biln(self):
        """Return a BILN string from the PolymerInfo generated at
        instantiation.
        
        :return: str"""
        biln = self.__to_biln()
        return biln

    # End of Converter pipeline functions.

############################################################
# End of converter.py
############################################################
