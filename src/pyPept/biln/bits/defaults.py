"""
A place where default BILN-related symbols and related information are defined.

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

# System libraries needed by this module.

# Third-party libraries needed by this module, e.g. numpy.

# Project-specific modules additionally needed.
from pyPept.biln.bits.abstractions import TitledBILNExample
from pyPept.monomerlib import MonomerConstants, MonomerLibrary

_installedMonoLib = MonomerLibrary()


# ----- Begin code for this module. -----

# Default delimiters of BILN
_defChainSeparator = "."
_defMonomerSeparator = "-"
_defNumAnnotationSeparator = ":"
_defNNAAStartDelim = "["
_defNNAAEndDelim = "]"
_defCompressCharacter = "x"


# Default threshold of linker/FA monomer fraction in a chain to be called HLE
_defLinkerFractionInHLE = 0.40

##########################################################################
# Functions and classes
##########################################################################

############################################################
class BILNConstants():
    """A placeholder for BI peptide line notation constants and defaults."""

    @staticmethod
    def GetDefaultMonomerSeparator():
        """
        Provides the expected/default monomer separator (character) in BILN.
        """
        return _defMonomerSeparator
    ########################################

    @staticmethod
    def GetDefaultChainSeparator():
        """
        Provides the expected/default chain separator (character) in BILN.

        >>> BILNConstants.GetDefaultChainSeparator()
        '.'
        """
        return _defChainSeparator
    ########################################

    @staticmethod
    def GetDefaultNumericAnnotationSeparator():
        """
        Provides the expected/default separator between chain numbers,
        monomer numbers, and monomer labels.

        >>> BILNConstants.GetDefaultNumericAnnotationSeparator()
        ':'
        """
        return _defNumAnnotationSeparator
    ########################################

    @staticmethod
    def GetDefaultNNAAStartDelimiter():
        """
        Provides the expected/default character that begins annotation of
        a non-natural amino acid or branch point.

        >>> BILNConstants.GetDefaultNNAAStartDelimiter()
        '['
        """
        return _defNNAAStartDelim
    ########################################

    @staticmethod
    def GetDefaultNNAAEndDelimiter():
        """
        Provides the expected/default character that ends annotation of
        a non-natural amino acid or branch point.

        >>> BILNConstants.GetDefaultNNAAEndDelimiter()
        ']'
        """
        return _defNNAAEndDelim
    ########################################

    @staticmethod
    def GetDefaultCompressedBILNDelimiter():
        """
        Provides the expected/default character that is used to indicate
        a consecutive non-breaking repeats of a monomer in a BILN.

        >>> BILNConstants.GetDefaultCompressedBILNDelimiter()
        'x'
        """
        return _defCompressCharacter
    ########################################

    @staticmethod
    def GetAcceptableMonomerCodes(monomer_library=None):
        """
        Provides a tuple of monomer IDs (strings) that are accepted as
        valid in BI Line Notation.

        :param monomer_library: a monomer library to be used as an override.
        :type monomer_library: MonomerLibrary

        >>> len(BILNConstants.GetAcceptableMonomerCodes())
        322
        """
        if monomer_library is not None:
            useLib = monomer_library
        else:
            useLib = _installedMonoLib  # default
        assert isinstance(useLib, MonomerLibrary)
        return tuple(useLib[MonomerConstants.attr_monomer_symbol])
    ########################################

    @staticmethod
    def GetMaximumRGroupNumber():
        """
        Provide the maximum R-group attachment point number tolerated by each
        monomer. Though four R-groups are theoretically possible, it is
        an observation that peptide monomers nearly always have R1- and R2-
        groups to form amide bonds, and a single R3 side-chain group.
        """
        return 3
    ########################################

    @staticmethod
    def GetBranchingRegex():
        """
        Provides a regular expression that can be used to detect loops or
        peptide chain-branching (e.g., HLE) attachment points."""
        return "\\(-?[1-9],-?[1-9]\\)"
        #return "\(-\d+,-?\d+)"  # Support arbitrary integer branchID!
    ########################################

    @staticmethod
    def GetFattyAcidRegex():
        """Provides a regular expression to detect fatty acids in a BILN."""
        return "C[012]?[0-9][D]?A"
    ########################################

    @staticmethod
    def _GetStandardAAs_1L(format="string"):
        _standard_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        if format == "string":
            return _standard_amino_acids
        # Otherwise as a tuple
        return tuple(l for l in _standard_amino_acids)
    ########################################

    # Define acceptable words used to define strategy for splitting mol to #chains
    labelSortLength = 'purelength'
    labelSortAACount = 'aacount'
    labelSortAAFrac = 'aaproportion'
    labelSortUnchanged = 'as-is'

    @staticmethod
    def GetDefaultChainSortCriteria():
        """Provides a sort criteria as defined in the BILNConstants.
        At present, this is to sort by chain length.
        """
        return BILNConstants.labelSortLength

    @staticmethod
    def GetInvalidSequenceExamples():
        """Provides a collection of INVALID BILN strings.

        :return: tuple of str

        >>> len(BILNConstants.GetInvalidSequenceExamples())
        14
        """
        return (
            TitledBILNExample(BILN="A--A-C",
                title="Non-library monomer example 1 (empty string)"),
            TitledBILNExample(BILN="A-MYMISTAKENMONOMER-C",
                title="Non-library monomer example 2"),
            TitledBILNExample(BILN="A-A-C.K",
                title="Ambiguous attachment example 1"),
            TitledBILNExample(BILN="H-Aib-E-G-L-V-R-G-R-K(1,3)-OH",
                title="Ambiguous attachment example 2"),
            TitledBILNExample(BILN="H-Y-R-Q-R-Y-NH2.C18DA-OEG(1,2)",
                title="Ambiguous attachment example 3"),
            TitledBILNExample(BILN="Y-K(1,3)-P-P-P-S-NH2." +
                "C18DA(2,2)-gGlu-OEG(1,2).Gly(2,1)",
                title="Non-library monomer example 3 ('Gly' should be 'G')"),
            TitledBILNExample(BILN="H-AC4C-Q",
                title="Non-library monomer example 4 (should be 'Ac4C')"),
            TitledBILNExample(BILN=\
                "C(1,3)-A-I-C(1,3)-L-V-NH2.bAcAPen(1,1)(2,2)",
                title="Valid pair branch but invalid bond indices."),
            TitledBILNExample(BILN="G(1,4)-K(1,3)-D(2,2).OEG2(1,2)",
                title="No OEG2 monomer, bond indices mismatch, no R4 for G"),

            # R-group assignment problems
            TitledBILNExample(BILN="Y-K(1,1)-G-Y-NH2.C20DA-gGlu-eLys(1,2)",
                title="R1 invalid when in middle of chain."),
            TitledBILNExample(BILN="Y-K(1,2)-G-Y-NH2.C20DA-gGlu-eLys(1,2)",
                title="R2 invalid when in middle of chain."),
            TitledBILNExample(BILN="K(1,2)-G-Y-NH2.C20DA-gGlu-eLys(1,2)",
                title="R2 invalid when at the front of a chain."),
            TitledBILNExample(BILN="Y-K(1,3)-G-Y-NH2.C20DA-gGlu-eLys(1,1)",
                title="R1 invalid when at the end of a chain."),
            TitledBILNExample(BILN="Y-C(1,3)-G-Y-C(1,1)",
                title="R1 invalid when at the end of a chain."),

        )
    ########################################

    @staticmethod
    def GetValidSequenceExamples():
        """Provides a collection of VALID BILN strings.

        :return: tuple of str

        >>> len(BILNConstants.GetValidSequenceExamples())
        5
        """
        return tuple(
            [TitledBILNExample(BILN=seq, title="Valid example %i" % i)
                for i, seq in enumerate(
            (
            "C(1,3)-A-A-A-C(1,3)",
            "C(1,1)-A-A-A-C(1,2)",
            "A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)-A-A-D(2,3)",
            "A-gGlu-gGlu-OEG-OEG-Y-gGlu-NMeL-Nle-Ac4C-P-R-S-NH2",
            "A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2)." + \
                "C-G-G-G(2,2).C-G-G-G(3,2)",
            ), start=1)])
    ########################################

    @staticmethod
    def GetAllSequenceExamples():
        """Combines the collections of invalid and valid BILN strings.

        :return: tuple of str

        >>> len(BILNConstants.GetAllSequenceExamples())
        19
        """
        return tuple(
            list(BILNConstants.GetInvalidSequenceExamples()) +
            list(BILNConstants.GetValidSequenceExamples())
            )
    ########################################


    @staticmethod
    def GetDefaultHLELinkerFraction():
        """Provides a default value for the fraction of monomers in a chain
        that must be present as FAs (fatty acids) or linker monomers.
        Should be positive value.

        >>> BILNConstants.GetDefaultHLELinkerFraction()
        0.4
        """
        return _defLinkerFractionInHLE
    ########################################

## end of BILNConstants class ##############################


## Verify that the module's interfaces work as the doctests demonstrate.
if __name__ == "__main__":

    import doctest, os, sys, __main__
    num_fail, numTests = doctest.testmod()
    modulePath = os.path.abspath(__main__.__file__)
    if num_fail > 0:
        print('%s : Expected functionality in doctests fail!' % modulePath)
        sys.exit(-1)
    elif numTests > 0:
        print("%s : All %i doctests passed." % (modulePath, numTests))
        sys.exit(0)
    else:
        print("%s : WARNING - No doctests defined." % modulePath)
        sys.exit(1)
