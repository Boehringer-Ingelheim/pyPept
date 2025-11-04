"""
A module dedicated to the validity of BILN strings and information extraction.

This module provides methods to confirm what a valid BILN string is, and
what sort of information can be extracted from a BILN.

Notes about BILN important to interpretation:
A monomer is expressed in the format <X>(a,b) , where:
    <X> is the monomer name
    a is an identifier of a cross-linked bond (inter- or intra-chain)
    b is the R-group where a bond occurs:
        b=1   means the N-terminus of a monomer
        b=2   means the C-terminus of a monomer
        b=3,4 means another R-group used for linkages

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
import logging   # For messages and debugging.
import re
from collections import defaultdict
import itertools

# Third-party libraries needed by this module, e.g. numpy.


# Project-specific modules additionally needed.
from pyPept.biln.bits.abstractions import (
    NumberedChain,
    NumberedMonomer, MonomerPair,
    InvalidMonomerName, AmbiguousBranchIDs,
    InvalidMonomerRGroup, InvalidChainBeginRGroup, InvalidChainEndRGroup,
    InvalidChainMiddleRGroup, RepeatRGroupMonomer, BILNMultiError,
    )
from pyPept.biln.bits.defaults import BILNConstants, _installedMonoLib
from pyPept.monomerlib import MonomerLibrary

# ----- Begin code for this module. -----


# Convenience function used in different parts of the parser.
sortMonomers = lambda m: (m.chain, m.number)    # Reliant on type-checking


##########################################################################
# Functions and classes
##########################################################################

############################################################
class BILNParser():
    """
    A class which provides information services and parsing on BILN strings.
    """

    ############################################################
    @staticmethod
    def IsStandardAA(seq):
        """Returns whether or not a given string is a monomer representing
        a standard amino acid. Here, the rule used is that standard AAs
        are single letters not including {BJOUXZ}.

        >>> BILNParser.IsStandardAA("A")
        True
        >>> BILNParser.IsStandardAA("K(1,3)")
        True
        >>> BILNParser.IsStandardAA("[K(1,3)]") # Occurs when reformatting.
        True
        >>> BILNParser.IsStandardAA("Lys(1,3)") # Violates one-letter rule.
        False
        """
        stripped = re.sub(BILNConstants.GetBranchingRegex(),
                    str(), seq)
        stripped = re.sub("\\" + BILNConstants.GetDefaultNNAAStartDelimiter(),
                    str(), stripped)
        stripped = re.sub("\\" + BILNConstants.GetDefaultNNAAEndDelimiter(),
                    str(), stripped)
        return stripped in BILNConstants._GetStandardAAs_1L()
    ############################################################

    ############################################################
    @staticmethod
    def GetPlainMonomerCode(seq):
        """Strips a monomer of potential branch and non-standard AA notation.

        :param seq: a lone monomer in its full BILN representation
        :type seq: str
        :return: str

        >>> BILNParser.GetPlainMonomerCode("[K(1,3)]")
        'K'
        >>> BILNParser.GetPlainMonomerCode("K(1,3)")
        'K'
        >>> BILNParser.GetPlainMonomerCode("[OEG(2,2)]")
        'OEG'
        >>> BILNParser.GetPlainMonomerCode("[OH]")
        'OH'
        >>> BILNParser.GetPlainMonomerCode("[C20DA]")
        'C20DA'
        >>> BILNParser.GetPlainMonomerCode("[Aib]")
        'Aib'
        """
        if seq is None:
            return None
        no_branch = re.sub(BILNConstants.GetBranchingRegex(),
                        str(), seq)
        no_branch = re.sub("\\" + BILNConstants.GetDefaultNNAAStartDelimiter(),
                        str(), no_branch)
        no_branch = re.sub("\\" + BILNConstants.GetDefaultNNAAEndDelimiter(),
                        str(), no_branch)
        return no_branch
    ############################################################

    ############################################################
    @staticmethod
    def IsValidMonomer(seq, raiseError=False, monomer_library=None,
        ):
        """For a given monomer string and its optional branching annotation,
            determine if the monomer name is valid and if the optional
            branching(s) and R-group attachment point indices are valid or not.

        >>> BILNParser.IsValidMonomer('A')
        True
        >>> BILNParser.IsValidMonomer('K')
        True
        >>> BILNParser.IsValidMonomer('K(1,1)')
        True
        >>> BILNParser.IsValidMonomer('K(1,2)')
        True
        >>> BILNParser.IsValidMonomer('K(1,3)')
        True
        >>> BILNParser.IsValidMonomer('K(1,4)')
        False
        >>> BILNParser.IsValidMonomer('K(3,3)')
        True
        >>> BILNParser.IsValidMonomer('K(4,-1)')
        False
        >>> BILNParser.IsValidMonomer('K(4,3)')
        True
        >>> BILNParser.IsValidMonomer('K(1,1)(2,2)')
        True
        >>> BILNParser.IsValidMonomer('K(1,1)(2,1)')
        False
        >>> BILNParser.IsValidMonomer('K(1,1)(2,2)(3,3)')
        True
        >>> BILNParser.IsValidMonomer('MYWRONGMONOMER(4,3)')
        False

        Checking of the types of errors that could be raised and caught.
        >>> BILNParser.IsValidMonomer('K(4,-1)', raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidMonomerRGroup: K(4,-1)
        >>> BILNParser.IsValidMonomer('K(1,1)(2,1)', raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.RepeatRGroupMonomer: K(1,1)(2,1)
        >>> BILNParser.IsValidMonomer('MYWRONGMONOMER(4,3)', raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidMonomerName: MYWRONGMONOMER(4,3)
        >>> BILNParser.IsValidMonomer('K(1,4)(2,5)', raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [InvalidMonomerRGroup('K(1,4)(2,5)'), InvalidMonomerRGroup('K(1,4)(2,5)')]
        """
        if monomer_library is not None:
            useLib = monomer_library
        else:
            useLib = _installedMonoLib  # default
        assert isinstance(useLib, MonomerLibrary)

        errors = list()
        validMonomers = BILNConstants.GetAcceptableMonomerCodes()
        basename = BILNParser.GetPlainMonomerCode(seq)
        if basename not in validMonomers:
            errors.append(InvalidMonomerName(seq))
        else: # monomer valid, any R-group error?
            branchTexts = re.findall(BILNConstants.GetBranchingRegex(), seq)
            if branchTexts: # R groups listed valid?
                brNums = [b.split(",")[-1].replace(")","") for b in branchTexts]
                if len(brNums) != len(set(brNums)):  # Duplicate R-group
                    errors.append(RepeatRGroupMonomer(seq))
                for brNum in brNums: # Check R-groups exist
                    if 'R'+brNum not in useLib.GetRGroups(basename):
                        errors.append(InvalidMonomerRGroup(seq))
        # Checked monomer name and R-groups, report valid or not
        if errors:
            if raiseError:
                if len(errors) > 1:
                    raise BILNMultiError(errors)
                else:
                    raise errors.pop()
            else:
                return False
        else:
            return True
    ############################################################

    ############################################################
    @staticmethod
    def GetInvalidMonomers(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            ):
        """Returns a list of monomers in the sequence that are invalid.

        >>> peps = BILNConstants.GetAllSequenceExamples()
        >>> [BILNParser.GetInvalidMonomers(p.BILN) for p in peps[:3]]
        [('',), ('MYMISTAKENMONOMER',), ()]
        """
        retMonomers = list()
        for mol in seq.split(chainSep):
            for mon in mol.split(monoSep):
                if not BILNParser.IsValidMonomer(mon):
                    retMonomers.append(mon)
        return tuple(retMonomers)
    ############################################################

    ############################################################
    @staticmethod
    def _GetBranchIDs(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            ):
        """
        Returns the branch annotation values (IDs) found in a given BILN.

        For example, if a chain contains only a single branch to another chain,
        one would expect a single value of 1 in a tuple. Additionally, if
        a peptide is a chained HLE (one HLE attached to a main chain, with a
        second HLE attached to the first HLE), one would expect the values of
        1 and 2 in a tuple.

        :return: a tuple of integer

        >>> for s in BILNConstants.GetValidSequenceExamples():
        ...     print(s.BILN[:80])
        ...     print(BILNParser._GetBranchIDs(s.BILN))
        C(1,3)-A-A-A-C(1,3)
        (1,)
        C(1,1)-A-A-A-C(1,2)
        (1,)
        A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)-A-A-D(2,3)
        (1, 2)
        A-gGlu-gGlu-OEG-OEG-Y-gGlu-NMeL-Nle-Ac4C-P-R-S-NH2
        ()
        A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2).C-G-G-G(2,2).C-G-G-G(3,2)
        (1, 2, 3)
        """
        branchPairs = re.findall(BILNConstants.GetBranchingRegex(), seq)
        branchNums = [int(pair[1]) for pair in branchPairs]
        return tuple(sorted(set(branchNums)))
    ############################################################

    ############################################################
    @staticmethod
    def SplitToNumberedChains(seq, start=1,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            orderCriteria=BILNConstants.GetDefaultChainSortCriteria(),
            logger=None):
        """
        Splits a BILN into sub-chains and assigns each one an ID number.

        :param start: What number to use for assigning chain IDs.
        :type start: int
        :param orderCriteria: How to order the chains. Options are
            defined in the BILNConstants.labelSortXX data members.
        :type orderCriteria: str
        :param logger: A logger that could provide messages (currently unused).
        :type logger: A logging.Logger instance
        :return: a list of biln.NumberedChain instances (named tuples)

        >>> peps = BILNConstants.GetValidSequenceExamples()
        >>> bilnSplits = list() # List of lists.
        >>> for pep in peps:
        ...     bilnSplits.append(BILNParser.SplitToNumberedChains(pep.BILN))
        >>> for counter, ex_chains in enumerate(zip(peps, bilnSplits), start=1):
        ...     if len(ex_chains[1]) > 2: # 3+ chains, show results.
        ...         print("Sequence %i %s" % (counter, ex_chains[0].BILN))
        ...         for chain in ex_chains[1]:
        ...             print(chain)
        ...         print('-----')
        Sequence 5 A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2).C-G-G-G(2,2).C-G-G-G(3,2)
        NumberedChain(chain=1, BILN='C-G-G-G-K(2,3)-K(3,3)-G-G(1,2)')
        NumberedChain(chain=2, BILN='A-S-D-F-K(1,3)-A-S-F')
        NumberedChain(chain=3, BILN='C-G-G-G(3,2)')
        NumberedChain(chain=4, BILN='C-G-G-G(2,2)')
        -----

        Different ordering schemes supported:
        >>> biln = "K-E-K(1,3).C16DA-Sar-Sar-K(1,3)-Sar-K-E-A"
        >>> BILNParser.SplitToNumberedChains(biln)  # default sort by #monomers
        [NumberedChain(chain=1, BILN='C16DA-Sar-Sar-K(1,3)-Sar-K-E-A'), NumberedChain(chain=2, BILN='K-E-K(1,3)')]
        >>> BILNParser.SplitToNumberedChains(biln,
        ...     orderCriteria=BILNConstants.labelSortLength) # explicit specification of default sort
        [NumberedChain(chain=1, BILN='C16DA-Sar-Sar-K(1,3)-Sar-K-E-A'), NumberedChain(chain=2, BILN='K-E-K(1,3)')]
        >>> BILNParser.SplitToNumberedChains(biln,
        ...     orderCriteria=BILNConstants.labelSortAACount)
        [NumberedChain(chain=1, BILN='C16DA-Sar-Sar-K(1,3)-Sar-K-E-A'), NumberedChain(chain=2, BILN='K-E-K(1,3)')]
        >>> BILNParser.SplitToNumberedChains(biln,
        ...     orderCriteria=BILNConstants.labelSortAAFrac)
        [NumberedChain(chain=1, BILN='K-E-K(1,3)'), NumberedChain(chain=2, BILN='C16DA-Sar-Sar-K(1,3)-Sar-K-E-A')]
        >>> BILNParser.SplitToNumberedChains(biln,
        ...     orderCriteria=BILNConstants.labelSortUnchanged) # do NOT change order
        [NumberedChain(chain=1, BILN='K-E-K(1,3)'), NumberedChain(chain=2, BILN='C16DA-Sar-Sar-K(1,3)-Sar-K-E-A')]
        >>> BILNParser.SplitToNumberedChains(
        ...     "OEG-OEG(2,2).OEG-K(2,3)-OEG(1,2).C-K(1,3)-D-E-D",
        ...     orderCriteria=BILNConstants.labelSortUnchanged) # do NOT change order
        [NumberedChain(chain=1, BILN='OEG-OEG(2,2)'), NumberedChain(chain=2, BILN='OEG-K(2,3)-OEG(1,2)'), NumberedChain(chain=3, BILN='C-K(1,3)-D-E-D')]
        """
        if seq == str(): return list()
        molChains = seq.split(chainSep)
        if orderCriteria == BILNConstants.labelSortLength:
            sizes_splitMols = sorted([(len(mol.split(monoSep)), mol)
                                for mol in molChains])
        elif orderCriteria == BILNConstants.labelSortAACount:
            sizes_splitMols = sorted([(
                len([res for res in mol.split(monoSep) if
                    BILNParser.IsStandardAA(res)]), mol)
                for mol in molChains])
        elif orderCriteria == BILNConstants.labelSortAAFrac:
            sizes_splitMols = sorted([(len([res for res in mol.split(monoSep)
                if BILNParser.IsStandardAA(res)]) / len(mol.split(monoSep)), mol)
                for mol in molChains])
        elif orderCriteria == BILNConstants.labelSortUnchanged:
            return [NumberedChain(i, ch)
                for (i,ch) in enumerate(molChains, start=1)]
        else:
            raise ValueError("Unaccepted ordering scheme: %s" % orderCriteria)

        largestMol = NumberedChain(start, sizes_splitMols.pop()[1]) # Back largest
        return [largestMol,] + \
                BILNParser.SplitToNumberedChains(
                chainSep.join([pair[1] for pair in sizes_splitMols]),
                start=start+1, chainSep=chainSep, monoSep=monoSep,
                orderCriteria=orderCriteria)
    ############################################################

    ############################################################
    @staticmethod
    def Length(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            logger=None,
            ):
        """Return the total number of residues in the given sequence.

        >>> for peptide in BILNConstants.GetValidSequenceExamples():
        ...     print(peptide.BILN)
        ...     print(BILNParser.Length(peptide.BILN))
        C(1,3)-A-A-A-C(1,3)
        5
        C(1,1)-A-A-A-C(1,2)
        5
        A-A-D(1,3)-A-A-K(2,3)-A-A.K(1,3)-A-A-D(2,3)
        12
        A-gGlu-gGlu-OEG-OEG-Y-gGlu-NMeL-Nle-Ac4C-P-R-S-NH2
        14
        A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2).C-G-G-G(2,2).C-G-G-G(3,2)
        24
        """
        chains = BILNParser.SplitToNumberedChains(
                    seq=seq, chainSep=chainSep, monoSep=monoSep)
        return sum([len(c.BILN.split(monoSep)) for c in chains])
    ############################################################


    ############################################################
    @staticmethod
    def HasValidBranchAnnotations(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            raiseError=False,
            ):
        """Reports if the branch annotations are valid for a given chain.
        Essentially, checks that all branch/join annotations come in pairs.

        >>> [BILNParser.HasValidBranchAnnotations(pep.BILN) for pep in
        ...     BILNConstants.GetInvalidSequenceExamples()]
        [True, True, True, False, False, True, True, False, False, False, False, False, False, False]
        >>> all([BILNParser.HasValidBranchAnnotations(pep.BILN) for pep in
        ...     BILNConstants.GetValidSequenceExamples()])
        True

        Detailed errors:
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[3].BILN
        >>> pep
        'H-Aib-E-G-L-V-R-G-R-K(1,3)-OH'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.AmbiguousBranchIDs: 1
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[4].BILN
        >>> pep
        'H-Y-R-Q-R-Y-NH2.C18DA-OEG(1,2)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.AmbiguousBranchIDs: 1
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[7].BILN
        >>> pep
        'C(1,3)-A-I-C(1,3)-L-V-NH2.bAcAPen(1,1)(2,2)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [AmbiguousBranchIDs(1), AmbiguousBranchIDs(2)]
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[8].BILN
        >>> pep
        'G(1,4)-K(1,3)-D(2,2).OEG2(1,2)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [AmbiguousBranchIDs(1), AmbiguousBranchIDs(2)]

        # Invalid R-group annotations during branching
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[-5].BILN
        >>> pep
        'Y-K(1,1)-G-Y-NH2.C20DA-gGlu-eLys(1,2)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidChainMiddleRGroup: Y-K(1,1)-G-Y-NH2
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[-4].BILN
        >>> pep
        'Y-K(1,2)-G-Y-NH2.C20DA-gGlu-eLys(1,2)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidChainMiddleRGroup: Y-K(1,2)-G-Y-NH2
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[-3].BILN
        >>> pep
        'K(1,2)-G-Y-NH2.C20DA-gGlu-eLys(1,2)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidChainBeginRGroup: K(1,2)-G-Y-NH2
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[-2].BILN
        >>> pep
        'Y-K(1,3)-G-Y-NH2.C20DA-gGlu-eLys(1,1)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidChainEndRGroup: C20DA-gGlu-eLys(1,1)
        >>> pep = BILNConstants.GetInvalidSequenceExamples()[-1].BILN
        >>> pep
        'Y-C(1,3)-G-Y-C(1,1)'
        >>> BILNParser.HasValidBranchAnnotations(pep, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidChainEndRGroup: Y-C(1,3)-G-Y-C(1,1)
        """
        brRegex = BILNConstants.GetBranchingRegex()
        errors = list()
        branchPairs = re.findall(brRegex, seq)
        branchNums = [int(pair[1]) for pair in branchPairs] # 2nd char in '(X,Y)'
        for branchID in sorted(set(branchNums)):
            if branchNums.count(branchID) != 2:
                errors.append(AmbiguousBranchIDs(branchID))
        # Now, check to be sure that branch annotation context is correct.
        #   Cannot have R-group R2 at front of chain, nor R1 at end of chain,
        #   nor R1/R2 in middle of chain.
        for chain in BILNParser.SplitToNumberedChains(seq,
                chainSep=chainSep, monoSep=monoSep):
            monomers = chain.BILN.split(monoSep)
            if len(monomers) == 1:  # Don't check special monomers such as
                                    # Cy7() or bAcAPen()()
                continue
            for pos, mon in enumerate(monomers, start=1):
                branches = re.findall(brRegex, mon)
                if branches:
                    rGrp = int(mon.split(",")[-1][0])
                    if pos == 1 and rGrp == 2 and chain.chain == 1:
                        errors.append(InvalidChainBeginRGroup(
                            "%s" % chain.BILN))
                    elif pos == len(monomers) and rGrp == 1:
                        errors.append(InvalidChainEndRGroup(
                            "%s" % chain.BILN))
                    elif len(monomers) > 2 and rGrp in (1,2) and \
                                                1 < pos < len(monomers):
                        errors.append(InvalidChainMiddleRGroup(
                            "%s" % chain.BILN))
        # End of error checking
        if errors:
            if raiseError:
                if len(errors) > 1:
                    raise BILNMultiError(errors)
                else:
                    raise errors.pop()
            else:
                return False
        else:
            return True
    ############################################################

    ########################################
    @staticmethod
    def GetNNAAMonomers(
            biln,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            ):
        """Retrieve the non-natural amino acid monomers in a BILN.

        :return: a tuple of NumberedMonomer objects.

        >>> for m in BILNParser.GetNNAAMonomers("4Pal-A-Aib-K-S-P"):
        ...     print(m)
        NumberedMonomer(chain=1, number=1, monomer='[4Pal]')
        NumberedMonomer(chain=1, number=3, monomer='[Aib]')
        """
        monomerObjs = BILNParser.GetNumberedBILN(seq=biln,
            chainSep=chainSep, monoSep=monoSep, returnFormat='monomerobj')
        return tuple([m for m in monomerObjs
            if not BILNParser.IsStandardAA(m.monomer)])
    ########################################

    ############################################################
    @staticmethod
    def GetNumberedBILN(seq, raiseError=False,
            start=1, resetMonomerNumbers=False,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            orderCriteria=BILNConstants.GetDefaultChainSortCriteria(),
            annoSep=BILNConstants.GetDefaultNumericAnnotationSeparator(),
            chainID=True,
            nonStdLeft=BILNConstants.GetDefaultNNAAStartDelimiter(),
            nonStdRight=BILNConstants.GetDefaultNNAAEndDelimiter(),
            monomerID=True,
            returnFormat="string",
            monomerContext=True,
            logger=None):
        """
        Annotate a raw BILN with chain and monomer number data.

        :param resetMonomerNumbers: Whether or not to start monomer numbering at the
                            given start value, for each chain.
        :param resetMonomerNumbers: bool
        :param returnFormat: How the number-annotated BILN should be returned.
                             Accepted values are 'string' or 'monomerobj'.
        :type returnFormat: str
        :param monomerContext: Whether or not monomers should be stripped to
                               just their name
        :type monomerContext: bool
        :return: Depends on return format value.
                 'string' provides a single string,
                 'monomerobj' provides an iterable of NumeredMonomer objects.

        >>> usePep = BILNConstants.GetValidSequenceExamples()[-1]
        >>> usePep.BILN
        'A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2).C-G-G-G(2,2).C-G-G-G(3,2)'

        Return in NumberedMonomer format:
        >>> retValue = BILNParser.GetNumberedBILN(usePep.BILN,
        ...                                         returnFormat='monomerobj')
        >>> type(retValue), len(retValue)
        (<class 'tuple'>, 24)
        >>> len(retValue) == BILNParser.Length(usePep.BILN)
        True

        Note the use of sort criteria for chain ordering, which is used in
        order to assist canonicalization of a large collection of BILNs.
        That means monomers may not come back in the order they were given
        in the original BILN:
        >>> [m for m in retValue if "[" in m.monomer][0]
        NumberedMonomer(chain=1, number=5, monomer='[K(2,3)]')

        By default, extra notation is given to non-standard peptide monomers,
        so as to easily see NNAAs or branch points. This can be turned off
        easily:
        >>> retValue = BILNParser.GetNumberedBILN(usePep.BILN,
        ...     returnFormat='monomerobj', nonStdLeft="", nonStdRight="")
        >>> [m for m in retValue if m.chain == 1 and m.number == 5][0]
        NumberedMonomer(chain=1, number=5, monomer='K(2,3)')

        Return in NumberedMonomer format, with context removed:
        >>> retValue = BILNParser.GetNumberedBILN(usePep.BILN,
        ...         returnFormat='monomerobj', monomerContext=False)
        >>> [m for m in retValue if m.chain == 1 and m.number == 5][0]
        NumberedMonomer(chain=1, number=5, monomer='K')

        Return as string, with options to turn on/off chains and monomer IDs:
        Again, note that the chain ordering can yield a result different from
        the original input BILN:
        >>> print(BILNParser.GetNumberedBILN(usePep.BILN,
        ...     chainID=False, monomerID=False))
        C-G-G-G-[K(2,3)]-[K(3,3)]-G-[G(1,2)].A-S-D-F-[K(1,3)]-A-S-F.C-G-G-[G(3,2)].C-G-G-[G(2,2)]
        >>> print(BILNParser.GetNumberedBILN(usePep.BILN)[:85])  # default behavior
        1:1:C-1:2:G-1:3:G-1:4:G-1:5:[K(2,3)]-1:6:[K(3,3)]-1:7:G-1:8:[G(1,2)].2:9:A-2:10:S-2:1
        """
        __acceptFmts = ("string", "monomerobj")
        if returnFormat not in __acceptFmts:
            raise ValueError(
                "returnFormat should be in %s" % " ".join(__acceptFmts))
        # This method is called during IsValidSequence(), it must self-resolve.
        retChains, resNum, asNumbered = list(), 0, list()
        chains = BILNParser.SplitToNumberedChains(seq=seq, start=start,
            chainSep=chainSep, monoSep=monoSep, orderCriteria=orderCriteria)
        for chain in chains:
            chainResidues = list()
            if resetMonomerNumbers and chainID:
                resNum = 0
            for mon in chain.BILN.split(monoSep):
                resNum += 1
                if len(mon) == 1:
                    useText = mon
                else:
                    if len(mon) >= 3 and \
                            mon[0] == nonStdLeft and mon[-1] == nonStdRight:
                        useText = mon
                    else: # Add standardizing notation
                        useText = nonStdLeft + mon + nonStdRight # Should be in standardize() or similar!

                if chainID:
                    chainResidues.append(annoSep.join(
                        (str(chain.chain), str(resNum), useText)))
                else:
                    if monomerID:
                        chainResidues.append(annoSep.join(
                            (str(resNum), useText)))
                    else: # Just getting flat, unified text
                        chainResidues.append(useText)
                # Store NumberedMonmer format in case requested.
                asNumbered.append(NumberedMonomer(chain=chain.chain,
                    number=resNum, monomer=useText))
            # Finished monomers in chain, store
            retChains.append(monoSep.join(chainResidues))
        # Finished re-constructing chains, now join all chains
        if returnFormat == 'monomerobj':
            if not monomerContext:
                asNumbered = [NumberedMonomer(chain=m.chain, number=m.number,
                    monomer=BILNParser.GetPlainMonomerCode(m.monomer))
                    for m in asNumbered]
            return tuple(asNumbered)
        else:
            return chainSep.join(retChains)
    ########################################



    ############################################################
    @staticmethod
    def GetCompressedMonomerDetectionRegex():
        """
        Provides a regular expression that can be used to detect consecutive
        non-branching repeats of a monomer in a BILN.
        """
        defDelim = BILNConstants.GetDefaultCompressedBILNDelimiter()
        if defDelim not in ("+", "*", "?"):
            return "\\d+" + defDelim
        else:
            return "\\d+\\" + defDelim
    ############################################################

    ############################################################
    @staticmethod
    def CompressBILN(seq,
            minRepeat=2,
            joinChar=BILNConstants.GetDefaultCompressedBILNDelimiter(),
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            nonStdLeft=BILNConstants.GetDefaultNNAAStartDelimiter(),
            nonStdRight=BILNConstants.GetDefaultNNAAEndDelimiter(),
            logger=None):
        """
        Compress a BILN by shortening the repeating non-branching monomer units
        down to a single unit along with its repeat length.
        Typical applications are to shorten HLE or linker strings, or to
        pass compressed BILNs on to other tools such as string searches in
        order to systematically identify all sequences with sequential repeats.

        :seealso: GetNumberedBILN

        >>> BILNParser.CompressBILN("A-A-A-A")
        '4xA'
        >>> BILNParser.CompressBILN("A-A-A-A-K(1,3)-A-A-A.OEG-OEG-OEG(1,2)")
        '4xA-K(1,3)-3xA.2xOEG-OEG(1,2)'

        Note the change in result if any NNAA or branch is part of the input:
        >>> BILNParser.CompressBILN("A-A-A-A-[K(1,3)]-A-A-A.OEG-OEG-OEG(1,2)")
        '4xA-[K(1,3)]-3xA.2x[OEG]-[OEG(1,2)]'

        >>> BILNParser.CompressBILN(minRepeat=3,
        ...     seq="A-A-A-A-K(1,3)-A-A-A.OEG-OEG-OEG(1,2)")
        '4xA-K(1,3)-3xA.OEG-OEG-OEG(1,2)'
        """
        passNSL = BILNConstants.GetDefaultNNAAStartDelimiter()
        passNSR = BILNConstants.GetDefaultNNAAEndDelimiter()
        if passNSL not in seq and passNSR not in seq:
            passNSL, passNSR = "", ""
        chains, retChains = list(), list()
        for chain in BILNParser.SplitToNumberedChains(
                seq=seq, chainSep=chainSep, monoSep=monoSep,
                orderCriteria=BILNConstants.labelSortUnchanged):

            chains.append([m.monomer for m in
                BILNParser.GetNumberedBILN(seq=chain.BILN,
                    chainSep=chainSep, monoSep=monoSep,
                    nonStdLeft=passNSL, nonStdRight=passNSR,
                    returnFormat="monomerobj")
                    ])
        for oneCh in chains:
            addChain = list()
            for symbol, replicates in itertools.groupby(oneCh):
                repCount = len([v for v in replicates])
                if repCount == 1:
                    addChain.append(symbol)
                else:
                    if repCount >= minRepeat:
                        addChain.append("%i%s%s" %
                                        (repCount, joinChar, symbol))
                    else: # Not enough replicates
                        addChain.append(monoSep.join([symbol]*repCount))
            retChains.append(addChain)
        retBILN = chainSep.join([monoSep.join(c) for c in retChains])
        return retBILN
    ############################################################

    ############################################################
    @staticmethod
    def DecompressBILN(seq,
            joinChar=BILNConstants.GetDefaultCompressedBILNDelimiter(),
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            logger=None):
        """
        De-compress a BILN by expanding repeating non-branching monomer units
        as indicated in the compressed BILN.
        Applications include quick BILN-to-molblock conversion, chemical
        registration systems, visual representation

        :param joinChar: What delimiter is used between the count and monomer.
        :type joinChar: str
        :return: str

        >>> BILNParser.DecompressBILN("8xA")
        'A-A-A-A-A-A-A-A'
        >>> BILNParser.DecompressBILN("A-K(1,3)-A.4xOEG-OEG(1,2)")
        'A-K(1,3)-A.OEG-OEG-OEG-OEG-OEG(1,2)'
        >>> compress = BILNParser.CompressBILN
        >>> decompress = BILNParser.DecompressBILN
        >>> usePeps = BILNConstants.GetValidSequenceExamples()
        >>> all(pep.BILN == decompress(compress(pep.BILN)) for pep in usePeps)
        True
        """
        chains = list()
        for chain in BILNParser.SplitToNumberedChains(
                seq=seq, chainSep=chainSep, monoSep=monoSep,
                orderCriteria=BILNConstants.labelSortUnchanged):
            monomers = list()
            for m in chain.BILN.split(monoSep):
                counter = re.findall("\\d+" + joinChar, m)
                if counter:
                    count = int(counter[0].replace(joinChar, ""))
                    monomers += count * [m.replace(counter[0], "")]
                else:
                    monomers.append(m)
            chains.append(monoSep.join(monomers))
        return chainSep.join(chains)
    ############################################################

    ############################################################
    @staticmethod
    def HasCompressedBILN(seq,
            joinChar=BILNConstants.GetDefaultCompressedBILNDelimiter(),
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            logger=None):
        """
        Indicate if a BILN contains a compressed representations.

        :param joinChar: What delimiter is used between the count and monomer.
        :type joinChar: str
        :return: str

        >>> BILNParser.HasCompressedBILN("8xA")
        True
        >>> BILNParser.HasCompressedBILN(BILNParser.DecompressBILN("8xA"))
        False

        >>> usePeps = BILNConstants.GetValidSequenceExamples() # No compressed BILNs
        >>> [BILNParser.HasCompressedBILN(pep.BILN) for pep in usePeps].count(True)
        0
        """
        for chain in BILNParser.SplitToNumberedChains(
                seq=seq, chainSep=chainSep, monoSep=monoSep,
                orderCriteria=BILNConstants.labelSortUnchanged):
            for m in chain.BILN.split(monoSep):
                counter = re.findall("\\d+" + joinChar, m)
                if counter:
                    return True
        return False
    ############################################################

    ############################################################
    @staticmethod
    def IsValidSequence(seq, raiseError=False,
            onlyMonoDictMonomers=True,
            joinChar=BILNConstants.GetDefaultCompressedBILNDelimiter(),
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            logger=None,
            ):
        """Returns whether a given BILN string is valid or not.

        :param seq: The BILN to test for validity.
        :type seq: str
        :param raiseException: whether or not to raise a ValueError when a
                               given BILN is invalid.
        :type raiseException: bool
        :param onlyMonoDictMonomers: Enforce monomer names are in dictionary.
        :type onlyMonoDictMonomers: bool
        :return: bool

        >>> any([BILNParser.IsValidSequence(pep.BILN)
        ...     for pep in BILNConstants.GetInvalidSequenceExamples()])
        False
        >>> all([BILNParser.IsValidSequence(pep.BILN)
        ...     for pep in BILNConstants.GetValidSequenceExamples()])
        False
        >>> all([BILNParser.IsValidSequence(pep.BILN, onlyMonoDictMonomers=False)
        ...     for pep in BILNConstants.GetValidSequenceExamples()])
        True

        Examples of all sequences with errors and the errors returned.
        There are 9 examples, currently.

        >>> invalids = BILNConstants.GetInvalidSequenceExamples()
        >>> invalids[1-1].BILN
        'A--A-C'
        >>> BILNParser.IsValidSequence(invalids[1-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidMonomerName: []
        >>> invalids[2-1].BILN
        'A-MYMISTAKENMONOMER-C'
        >>> BILNParser.IsValidSequence(invalids[2-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidMonomerName: [MYMISTAKENMONOMER]
        >>> invalids[3-1].BILN
        'A-A-C.K'
        >>> BILNParser.IsValidSequence(invalids[3-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.AmbiguousBranchIDs: A-A-C.K

        Currently, OH and NH2 capping monomers are not provided in the
        default library:
        >>> invalids[4-1].BILN
        'H-Aib-E-G-L-V-R-G-R-K(1,3)-OH'
        >>> BILNParser.IsValidSequence(invalids[4-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [AmbiguousBranchIDs(1), InvalidMonomerName('[OH]')]

        >>> invalids[6-1].BILN
        'Y-K(1,3)-P-P-P-S-NH2.C18DA(2,2)-gGlu-OEG(1,2).Gly(2,1)'
        >>> BILNParser.IsValidSequence(invalids[6-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [InvalidMonomerName('[NH2]'), InvalidMonomerName('[C18DA(2,2)]'), InvalidMonomerName('[OEG(1,2)]'), InvalidMonomerName('[Gly(2,1)]')]

        >>> invalids[7-1].BILN
        'H-AC4C-Q'
        >>> BILNParser.IsValidSequence(invalids[7-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.InvalidMonomerName: [AC4C]

        >>> invalids[8-1].BILN
        'C(1,3)-A-I-C(1,3)-L-V-NH2.bAcAPen(1,1)(2,2)'
        >>> BILNParser.IsValidSequence(invalids[8-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [BILNMultiError([AmbiguousBranchIDs(1), AmbiguousBranchIDs(2)]), InvalidMonomerName('[NH2]'), InvalidMonomerName('[bAcAPen(1,1)(2,2)]')]

        >>> invalids[9-1].BILN
        'G(1,4)-K(1,3)-D(2,2).OEG2(1,2)'
        >>> BILNParser.IsValidSequence(invalids[9-1].BILN, raiseError=True)
        Traceback (most recent call last):
        ...
        pyPept.biln.bits.abstractions.BILNMultiError: [BILNMultiError([AmbiguousBranchIDs(1), AmbiguousBranchIDs(2)]), InvalidMonomerRGroup('[G(1,4)]'), InvalidMonomerName('[OEG2(1,2)]')]
        """
        if BILNParser.HasCompressedBILN(
                seq=seq, joinChar=joinChar, chainSep=chainSep,
                monoSep=monoSep, logger=logger):
            return BILNParser.IsValidSequence(seq=BILNParser.DecompressBILN(
                    seq=seq, joinChar=joinChar, chainSep=chainSep,
                    monoSep=monoSep, logger=logger),
                joinChar=joinChar, chainSep=chainSep,
                monoSep=monoSep, logger=logger)

        # Step 1 -- valid branch notations?
        errors = list()
        try:
            BILNParser.HasValidBranchAnnotations(
                seq=seq, chainSep=chainSep, monoSep=monoSep, raiseError=True)
        except Exception as e:
            errors.append(e)
        # Step 2 -- monomer checking
        if onlyMonoDictMonomers:
            asMonomers = BILNParser.GetNumberedBILN(seq=seq,
                chainSep=chainSep, monoSep=monoSep, returnFormat="monomerobj")
            for m in asMonomers:
                try:
                    BILNParser.IsValidMonomer(m.monomer, raiseError=True)
                except Exception as e:
                    errors.append(e)
        # Step 3 -- does branch IDs detected match number of branches in mol?
        #   Needed to catch something like: A-A-C.K
        if len(BILNParser._GetBranchIDs(seq)) == 0 and \
            len(BILNParser.SplitToNumberedChains(seq)) > 1:
                errors.append(AmbiguousBranchIDs(seq))
        # Now final error handling.
        if errors:
            if raiseError:
                if len(errors) > 1:
                    raise BILNMultiError(errors)
                else:
                    raise errors.pop()
            else:
                if isinstance(logger, logging.Logger):
                    if len(errors) > 1:
                        logger.warning(BILNMultiError(errors))
                    else:
                        logger.warning(errors.pop())
                return False
        else:
            return True
    ## end of IsValidSequence static method ####################

    ############################################################
    @staticmethod
    def GetMonomerBranchIDs(seq):
        """
        Obtains the branch ID(s) from a branch-inclusive monomer.

        :param seq: a lone monomer in full BILN representation.
        :type seq: str
        :return: int or tuple(int) or None

        >>> BILNParser.GetMonomerBranchIDs("[K(1,3)]")
        1
        >>> BILNParser.GetMonomerBranchIDs("K(1,3)")
        1
        >>> BILNParser.GetMonomerBranchIDs("[OEG(2,2)]")
        2
        >>> BILNParser.GetMonomerBranchIDs("OEG") == None
        True
        >>> BILNParser.GetMonomerBranchIDs("bAcAPen(1,1)(2,2)")
        (1, 2)
        """
        retValue = None
        matches = re.findall(BILNConstants.GetBranchingRegex(), seq)
        if len(matches) == 1:
            retValue = int(matches.pop()[1]) # Assuming branch ID < 10
        elif len(matches) > 1:
            retValue = tuple([int(m[1]) for m in matches])
        return retValue
    ############################################################

    ############################################################
    @staticmethod
    def GetMonomerBranchRGroup(seq):
        """
        Obtains the R-group used in a branch-inclusive monomer.

        :param seq: a lone monomer in full BILN representation.
        :type seq: str
        :return: int or tuple(int) or None

        >>> BILNParser.GetMonomerBranchRGroup("[K(1,3)]")
        3
        >>> BILNParser.GetMonomerBranchRGroup("K(1,3)")
        3
        >>> BILNParser.GetMonomerBranchRGroup("[OEG(2,2)]")
        2
        >>> BILNParser.GetMonomerBranchRGroup("OEG") == None
        True
        >>> BILNParser.GetMonomerBranchRGroup("bAcAPen(1,1)(2,2)")
        (1, 2)
        """
        retValue = None
        matches = re.findall(BILNConstants.GetBranchingRegex(), seq)
        if len(matches) == 1:
            retValue = int(matches.pop()[3]) # Assuming branch ID < 10
        elif len(matches) > 1:
            retValue = tuple([int(m[3]) for m in matches])
        return retValue
    ############################################################

    ############################################################
    @staticmethod
    def GetMainPeptide(seq, minLength=3,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            monomerContext=True,
            orderCriteria=BILNConstants.GetDefaultChainSortCriteria(),
            logger=None,
            ):
        """Returns what is interpreted to be the main peptide chain
        of a molecule.

        :param seq: The BILN from which to extract the main peptide.
        :type seq: str
        :param minLength: The minimum length this peptide must be. A
                          heuristic of 3 is currently used.
        :param orderCriteria: How to order the chains. Options are
            defined in the BILNConstants.labelSortXX data members.
        :type orderCriteria: str
        :return: a NumberedChain representing the best guess of a main peptide.

        >>> validPeptides = BILNConstants.GetValidSequenceExamples()
        >>> for peptide in validPeptides:
        ...     print(BILNParser.GetMainPeptide(peptide.BILN))
        NumberedChain(chain=1, BILN='C(1,3)-A-A-A-C(1,3)')
        NumberedChain(chain=1, BILN='C(1,1)-A-A-A-C(1,2)')
        NumberedChain(chain=1, BILN='A-A-D(1,3)-A-A-K(2,3)-A-A')
        NumberedChain(chain=1, BILN='A-gGlu-gGlu-OEG-OEG-Y-gGlu-NMeL-Nle-Ac4C-P-R-S-NH2')
        NumberedChain(chain=1, BILN='C-G-G-G-K(2,3)-K(3,3)-G-G(1,2)')

        Strip chain context, such as when used in sequence alignment:
        >>> validPeptides[-1].BILN
        'A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2).C-G-G-G(2,2).C-G-G-G(3,2)'
        >>> BILNParser.GetMainPeptide(validPeptides[-1].BILN, monomerContext=False)
        NumberedChain(chain=1, BILN='C-G-G-G-K-K-G-G')

        Changing criteria in order to perceive main chain differently:
        >>> biln = "K-E-K(1,3).C16DA-Sar-Sar-K(1,3)-Sar-K-E-A"
        >>> BILNParser.GetMainPeptide(biln)
        NumberedChain(chain=1, BILN='C16DA-Sar-Sar-K(1,3)-Sar-K-E-A')
        >>> BILNParser.GetMainPeptide(biln, orderCriteria=BILNConstants.labelSortAAFrac)
        NumberedChain(chain=1, BILN='K-E-K(1,3)')
        """
        # Some state recording variables
        longest, longestNatAAs, largestFrac, doUpdate = None, None, None, False

        chains = BILNParser.SplitToNumberedChains(seq=seq, chainSep=chainSep,
            monoSep=monoSep, orderCriteria=orderCriteria)
        for chain in chains:
            doUpdate, chainLength = False, BILNParser.Length(chain.BILN)
            # Otherwise, already a longest chain, compare
            numNatAAs = len([m for m in chain.BILN.split(monoSep)
                            if BILNParser.IsStandardAA(m)])
            fracAAs = numNatAAs / chainLength

            if chainLength >= minLength: # Consider it as potential longest
                if longest is None:
                    doUpdate = True
                else:
                    if orderCriteria == BILNConstants.labelSortAAFrac:
                        if fracAAs > largestFrac:
                            doUpdate = True
                    elif numNatAAs > longestNatAAs:
                        doUpdate = True

                if doUpdate:
                    longest = chain
                    longestNatAAs, largestFrac = numNatAAs, fracAAs
        if monomerContext is False:
            return NumberedChain(chain=longest.chain,
                BILN=monoSep.join(
                    [BILNParser.GetPlainMonomerCode(m) for m in
                        longest.BILN.split(monoSep)]))
        else:
            return longest
    ############################################################

    ############################################################
    @staticmethod
    def ReconstructSequence(numberedMonomers,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            logger=None):
        """Rebuild a BILN from the NumberedMonomer objects given.

        :param numberedMonomers: a collection of NumberedMonomer objects.
        :type numberedMonomers: list, tuple, or set
        :return: str

        >>> from pyPept.biln import BILNConstants, BILNParser
        >>> get = BILNParser.GetNumberedBILN
        >>> recon = BILNParser.ReconstructSequence
        >>> all([get(example.BILN, chainID=False, monomerID=False) ==
        ...      recon(get(example.BILN, returnFormat="monomerobj"))
        ...      for example in BILNConstants.GetValidSequenceExamples()])
        True

        It is possible to reconstruct a BILN for fragments, that is, in which
        the sequence is not fully well-formed:
        >>> biln = "P-E-P-T-K(1,3)-D-E.OEG-OEG-K(2,3)-OEG-OEG(1,2).OEG(2,2)"
        >>> BILNParser.IsValidSequence(biln, onlyMonoDictMonomers=False)
        True
        >>> mons = BILNParser.GetNumberedBILN(biln, returnFormat='monomerobj',
        ...     nonStdLeft="", nonStdRight="")
        >>> len(mons)
        13
        >>> BILNParser.ReconstructSequence([m for m in mons if m.chain == 1])
        'P-E-P-T-K(1,3)-D-E'
        >>> BILNParser.ReconstructSequence([m for m in mons if m.chain != 1])
        'OEG-OEG-K(2,3)-OEG-OEG(1,2).OEG(2,2)'
        """
        chainBILNs, chainToResidues = list(), defaultdict(list)
        for res in numberedMonomers:
            assert isinstance(res, NumberedMonomer), \
                "Got %s (input was %s)" % (type(res), str(numberedMonomers))
            chainToResidues[res.chain].append(res)
        for chain in sorted(list(chainToResidues.keys())):
            chainResidues = [res.monomer for res in
                sorted(chainToResidues[chain], key=lambda r: r.number)]
            chainBILNs.append(monoSep.join(chainResidues))
        return chainSep.join(chainBILNs)
    ############################################################

    ############################################################
    @staticmethod
    def StripContext(
            biln,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            ):
        """Strip the branching/looping context from a given BILN.

        >>> BILNParser.StripContext("P-E-P-T-I-K(1,3)-E.OEG-OEG(1,2)")
        'P-E-P-T-I-K-E.OEG-OEG'
        >>> BILNParser.StripContext(
        ...     "P-E(1,3)-P-T-I-K(1,3)-K(2,3)-E.OEG-OEG(2,2)")
        'P-E-P-T-I-K-K-E.OEG-OEG'
        """
        numberedMonomers = BILNParser.GetNumberedBILN(biln,
            chainSep=chainSep, monoSep=monoSep, returnFormat="monomerobj",
            monomerContext=False)
        return BILNParser.ReconstructSequence(numberedMonomers,
            chainSep=chainSep, monoSep=monoSep)
    ############################################################

    ############################################################
    @staticmethod
    def _BILNToJoinPointNumberedMonomers(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            orderCriteria=BILNConstants.GetDefaultChainSortCriteria(),
            resetMonomerNumbers=False,
            nonStdLeft=BILNConstants.GetDefaultNNAAStartDelimiter(),
            nonStdRight=BILNConstants.GetDefaultNNAAEndDelimiter(),
            logger=None):
        """
        Internal-use method to retrieve annotated monomers from a BILN.
        Used by GetBranchMonomers() and GetIntraChainCycleMonomers().

        Potential liability:
        This method does *NOT* check for validity of branch annotations or
            validity of monomers (e.g., names and R-groups).

        :return: a defaultdict(list) of branch ID to monomer 2-tuples.

        >>> joinPoints = BILNParser._BILNToJoinPointNumberedMonomers(
        ...     seq="A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2)." +
        ...          "C-G-G-G(2,2).C-G-G-G(3,2)")
        >>> import pprint
        >>> pprint.pprint(joinPoints)
        defaultdict(<class 'list'>,
                    {1: [NumberedMonomer(chain=1, number=8, monomer='[G(1,2)]'),
                         NumberedMonomer(chain=2, number=13, monomer='[K(1,3)]')],
                     2: [NumberedMonomer(chain=1, number=5, monomer='[K(2,3)]'),
                         NumberedMonomer(chain=4, number=24, monomer='[G(2,2)]')],
                     3: [NumberedMonomer(chain=1, number=6, monomer='[K(3,3)]'),
                         NumberedMonomer(chain=3, number=20, monomer='[G(3,2)]')]})
        >>> joinPoints = BILNParser._BILNToJoinPointNumberedMonomers(
        ...     seq="C(1,3)-NMeI-V-T-P-I-C(2,3)-L-V-NH2.bAcAPen(1,1)(2,2)")
        >>> pprint.pprint(joinPoints)
        defaultdict(<class 'list'>,
                    {1: [NumberedMonomer(chain=1, number=1, monomer='[C(1,3)]'),
                         NumberedMonomer(chain=2, number=11, monomer='[bAcAPen(1,1)]')],
                     2: [NumberedMonomer(chain=1, number=7, monomer='[C(2,3)]'),
                         NumberedMonomer(chain=2, number=11, monomer='[bAcAPen(2,2)]')]})
        """
        monomers = BILNParser.GetNumberedBILN(seq=seq,
            resetMonomerNumbers=resetMonomerNumbers,
            chainSep=chainSep, monoSep=monoSep, returnFormat='monomerobj',
            orderCriteria=orderCriteria)
        monomers = [m for m in monomers if
            re.findall(BILNConstants.GetBranchingRegex(), m.monomer)]
        # Separate monomers into the branches they represent.
        branchID_to_monomers = defaultdict(list)
        while monomers:
            m = monomers.pop()
            branchIDs = BILNParser._GetBranchIDs(m.monomer)
            monBasename = BILNParser.GetPlainMonomerCode(m.monomer)
            for brID in branchIDs:
                newM = re.findall("\\("+str(brID)+",[1-9]\\)", m.monomer).pop()
                newM = nonStdLeft  + monBasename + newM + nonStdRight
                branchID_to_monomers[brID].append(NumberedMonomer(
                    chain=m.chain, number=m.number, monomer=newM,))
        # Sanity check
        failedIDs = [ID for ID,monomers in branchID_to_monomers.items() \
                        if len(monomers) != 2]
        if failedIDs:
            raise ValueError(("Unexpected monomer counts for %s: " % seq) +
                ";".join(["branch %i has %i monomer(s)" % (ID,
                len(branchID_to_monomers[ID])) for ID in failedIDs]) + ".")
        for brID, monomers in branchID_to_monomers.items():
            branchID_to_monomers[brID] = sorted(monomers)
        return branchID_to_monomers
    ############################################################

    ############################################################
    @staticmethod
    def GetBranchMonomers(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            orderCriteria=BILNConstants.GetDefaultChainSortCriteria(),
            monomerContext=True,
            resetMonomerNumbers=False,
            nonStdLeft=BILNConstants.GetDefaultNNAAStartDelimiter(),
            nonStdRight=BILNConstants.GetDefaultNNAAEndDelimiter(),
            logger=None):
        """
        Obtains pairs of monomers that form inter-chain branches.
        Does not handle any type of intra-chain branching/looping.

        :param seq: The full peptide molecule BILN to analyze for branches.
        :type seq: str
        :param monomerContext: If the context of a branching monomer should
                               be retained or not.
        :type monomerContext: bool
        :param resetMonomerNumbers: if monomers in each chain should be
                                    numbered from 1.
        :type resetMonomerNumbers: bool
        :return: a tuple of MonomerPair objects.

        When multiple branches are present, they are sorted by:
        1. branch ID,
        2. monomer chain ID,
        3. intra-chain monomer number.

        >>> peps = BILNConstants.GetValidSequenceExamples()
        >>> for pep in peps:
        ...     if len(BILNParser.GetBranchMonomers(pep.BILN)) > 0:
        ...         print(pep.BILN.replace("-",""))
        ...         print(len(BILNParser.GetBranchMonomers(pep.BILN)))
        ...         for mpair in BILNParser.GetBranchMonomers(pep.BILN):
        ...             print("Branch pair:")
        ...             print(mpair.monomer1)
        ...             print(mpair.monomer2)
        ...         print("-"*40)
        AAD(1,3)AAK(2,3)AA.K(1,3)AAD(2,3)
        2
        Branch pair:
        NumberedMonomer(chain=1, number=3, monomer='[D(1,3)]')
        NumberedMonomer(chain=2, number=9, monomer='[K(1,3)]')
        Branch pair:
        NumberedMonomer(chain=1, number=6, monomer='[K(2,3)]')
        NumberedMonomer(chain=2, number=12, monomer='[D(2,3)]')
        ----------------------------------------
        ASDFK(1,3)ASF.CGGGK(2,3)K(3,3)GG(1,2).CGGG(2,2).CGGG(3,2)
        3
        Branch pair:
        NumberedMonomer(chain=1, number=8, monomer='[G(1,2)]')
        NumberedMonomer(chain=2, number=13, monomer='[K(1,3)]')
        Branch pair:
        NumberedMonomer(chain=1, number=5, monomer='[K(2,3)]')
        NumberedMonomer(chain=4, number=24, monomer='[G(2,2)]')
        Branch pair:
        NumberedMonomer(chain=1, number=6, monomer='[K(3,3)]')
        NumberedMonomer(chain=3, number=20, monomer='[G(3,2)]')
        ----------------------------------------

        >>> print(BILNConstants.GetValidSequenceExamples()[5-1].BILN)
        A-S-D-F-K(1,3)-A-S-F.C-G-G-G-K(2,3)-K(3,3)-G-G(1,2).C-G-G-G(2,2).C-G-G-G(3,2)

        Option to drop a monomer context:
        >>> branches = BILNParser.GetBranchMonomers(
        ...     BILNConstants.GetValidSequenceExamples()[5-1].BILN,
        ...     monomerContext=False)
        >>> for mpair in branches:
        ...     print("Branch pair:")
        ...     print(mpair.monomer1)
        ...     print(mpair.monomer2)
        Branch pair:
        NumberedMonomer(chain=1, number=8, monomer='G')
        NumberedMonomer(chain=2, number=13, monomer='K')
        Branch pair:
        NumberedMonomer(chain=1, number=5, monomer='K')
        NumberedMonomer(chain=4, number=24, monomer='G')
        Branch pair:
        NumberedMonomer(chain=1, number=6, monomer='K')
        NumberedMonomer(chain=3, number=20, monomer='G')

        Option to reset monomer numbers in each chain:
        >>> branches = BILNParser.GetBranchMonomers(
        ...     BILNConstants.GetValidSequenceExamples()[5-1].BILN,
        ...     monomerContext=False, resetMonomerNumbers=True)
        >>> for mpair in branches:
        ...     print("Branch pair:")
        ...     print(mpair.monomer1)
        ...     print(mpair.monomer2)
        Branch pair:
        NumberedMonomer(chain=1, number=8, monomer='G')
        NumberedMonomer(chain=2, number=5, monomer='K')
        Branch pair:
        NumberedMonomer(chain=1, number=5, monomer='K')
        NumberedMonomer(chain=4, number=4, monomer='G')
        Branch pair:
        NumberedMonomer(chain=1, number=6, monomer='K')
        NumberedMonomer(chain=3, number=4, monomer='G')
        """
        branchID_to_monomers = BILNParser._BILNToJoinPointNumberedMonomers(
            seq=seq, chainSep=chainSep, monoSep=monoSep,
            orderCriteria=orderCriteria,
            resetMonomerNumbers=resetMonomerNumbers,
            nonStdLeft=nonStdLeft, nonStdRight=nonStdRight,
            logger=logger)
        # Known dict values are 2-tuples, keep only inter-branch
        dropIDs = list()
        for ID, monomers in branchID_to_monomers.items():
            chains = (monomers[0].chain, monomers[1].chain)
            if len(set(chains)) == 1:
                dropIDs.append(ID)
        for ID in dropIDs:
            del branchID_to_monomers[ID]
        if monomerContext is False:  # drop annotation from monomers
            for branchID in branchID_to_monomers.keys():
                branchID_to_monomers[branchID] = [
                    NumberedMonomer(chain=m.chain, number=m.number,
                        monomer=BILNParser.GetPlainMonomerCode(m.monomer))
                        for m in branchID_to_monomers[branchID]]
        # Sort to ensure lowest chain first, more natural interpretation.
        for branchID, monomers in branchID_to_monomers.items():
            branchID_to_monomers[branchID] = sorted(monomers, key=sortMonomers)
        return tuple([MonomerPair(
            monomer1=branchID_to_monomers[ID][0],
            monomer2=branchID_to_monomers[ID][1]) for
            ID in sorted(branchID_to_monomers.keys())])
    ############################################################

    @staticmethod
    def GetIntraChainCycleMonomers(seq,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            orderCriteria=BILNConstants.GetDefaultChainSortCriteria(),
            monomerContext=True,
            resetMonomerNumbers=False,
            nonStdLeft=BILNConstants.GetDefaultNNAAStartDelimiter(),
            nonStdRight=BILNConstants.GetDefaultNNAAEndDelimiter(),
            logger=None):
        """
        Obtains pairs of monomers that form intra-chain loops.
        Does not handle any type of inter-chain branching/looping.

        :see: GetBranchMonomers()
        :param seq: The full peptide molecule BILN to analyze for loops.
        :type seq: str
        :return: a tuple of MonomerPair objects.

        >>> peps = BILNConstants.GetValidSequenceExamples()
        >>> for pep in peps:
        ...     if len(BILNParser.GetIntraChainCycleMonomers(pep.BILN)) > 0:
        ...         print(pep.BILN.replace("-",""))
        ...         print(len(BILNParser.GetIntraChainCycleMonomers(pep.BILN)))
        ...         for mpair in BILNParser.GetIntraChainCycleMonomers(pep.BILN):
        ...             print(mpair.monomer1)
        ...             print(mpair.monomer2)
        C(1,3)AAAC(1,3)
        1
        NumberedMonomer(chain=1, number=1, monomer='[C(1,3)]')
        NumberedMonomer(chain=1, number=5, monomer='[C(1,3)]')
        C(1,1)AAAC(1,2)
        1
        NumberedMonomer(chain=1, number=1, monomer='[C(1,1)]')
        NumberedMonomer(chain=1, number=5, monomer='[C(1,2)]')

        Example of dropping monomer context:
        >>> peps = BILNConstants.GetValidSequenceExamples()
        >>> for pep in peps:
        ...     if len(BILNParser.GetIntraChainCycleMonomers(pep.BILN)) > 0:
        ...         print(pep.BILN.replace("-",""))
        ...         print(len(BILNParser.GetIntraChainCycleMonomers(pep.BILN)))
        ...         for mpair in BILNParser.GetIntraChainCycleMonomers(
        ...                 pep.BILN, monomerContext=False):
        ...             print(mpair.monomer1)
        ...             print(mpair.monomer2)
        ...         break
        C(1,3)AAAC(1,3)
        1
        NumberedMonomer(chain=1, number=1, monomer='C')
        NumberedMonomer(chain=1, number=5, monomer='C')

        Example of monomer numbering in second and later chains:
        >>> exampleBILN = "D-F-A-S-D-F-A-S-K(1,3)-F.C-E(2,3)-C-E(2,3)-OEG(1,2)"
        >>> for mpair in BILNParser.GetIntraChainCycleMonomers(exampleBILN,
        ...         resetMonomerNumbers=False):
        ...     print(mpair.monomer1)
        ...     print(mpair.monomer2)
        NumberedMonomer(chain=2, number=12, monomer='[E(2,3)]')
        NumberedMonomer(chain=2, number=14, monomer='[E(2,3)]')

        >>> for mpair in BILNParser.GetIntraChainCycleMonomers(exampleBILN,
        ...         resetMonomerNumbers=True):
        ...     print(mpair.monomer1)
        ...     print(mpair.monomer2)
        NumberedMonomer(chain=2, number=2, monomer='[E(2,3)]')
        NumberedMonomer(chain=2, number=4, monomer='[E(2,3)]')

        Impact of chain ordering on result:
        >>> biln = "P-E-P-K(1,3)-I.C20DA-C(2,3)-OEG-OEG-OEG-OEG-C(2,3)-OEG(1,2)"
        >>> for monopair in BILNParser.GetIntraChainCycleMonomers(biln):
        ...     print(monopair)
        MonomerPair(monomer1=NumberedMonomer(chain=1, number=2, monomer='[C(2,3)]'), monomer2=NumberedMonomer(chain=1, number=7, monomer='[C(2,3)]'))

        >>> for monopair in BILNParser.GetIntraChainCycleMonomers(biln,
        ...         orderCriteria=BILNConstants.labelSortAAFrac):
        ...     print(monopair)
        MonomerPair(monomer1=NumberedMonomer(chain=2, number=7, monomer='[C(2,3)]'), monomer2=NumberedMonomer(chain=2, number=12, monomer='[C(2,3)]'))
        """
        branchID_to_monomers = BILNParser._BILNToJoinPointNumberedMonomers(
            seq=seq, chainSep=chainSep, monoSep=monoSep,
            orderCriteria=orderCriteria,
            resetMonomerNumbers=resetMonomerNumbers,
            nonStdLeft=nonStdLeft, nonStdRight=nonStdRight,
            logger=logger)
        # Known dict values are 2-tuples, keep only intra-branch
        dropIDs = list()
        for ID, monomers in branchID_to_monomers.items():
            chains = (monomers[0].chain, monomers[1].chain)
            if len(set(chains)) != 1:
                dropIDs.append(ID)
        for ID in dropIDs:
            del branchID_to_monomers[ID]
        if monomerContext is False:  # drop annotation from monomers
            for branchID in branchID_to_monomers.keys():
                branchID_to_monomers[branchID] = [
                    NumberedMonomer(chain=m.chain, number=m.number,
                        monomer=BILNParser.GetPlainMonomerCode(m.monomer))
                        for m in branchID_to_monomers[branchID]]
        # Sort to ensure lowest chain first, more natural interpretation.
        for branchID, monomers in branchID_to_monomers.items():
            branchID_to_monomers[branchID] = sorted(monomers, key=sortMonomers)
        return tuple([MonomerPair(monomer1=pair[0], monomer2=pair[1]) for
            pair in branchID_to_monomers.values()])
    ############################################################

    ############################################################
    @staticmethod
    def GetIntraChainPath(numberedChain, monomer1, monomer2,
            chainSep=BILNConstants.GetDefaultChainSeparator(),
            monoSep=BILNConstants.GetDefaultMonomerSeparator(),
            keepEnds=True,
            logger=None):
        """
        Obtain the path of monomers between the two given monomers.

        :param numberedChain: the chain to look through for monomers.
        :type numberedChain: NumberedChain, with same chain ID as monomers
        :param monomer1: the first monomer to use in a path.
        :type monomer1: NumberedMonomer
        :param keepEnds: whether or not to keep the given monomers in the path
        :type keepEnds: bool
        :return: a tuple of NumberedMonomer objects.

        >>> usePep = "D-F-A-S-D-F-A-S-K(1,3)-F.C-E(2,3)-C-E(2,3)-OEG(1,2)"

        >>> chains = BILNParser.SplitToNumberedChains(usePep)
        >>> chains
        [NumberedChain(chain=1, BILN='D-F-A-S-D-F-A-S-K(1,3)-F'), NumberedChain(chain=2, BILN='C-E(2,3)-C-E(2,3)-OEG(1,2)')]

        Be careful with monomer numbering numbering:
        >>> mPairs = BILNParser.GetIntraChainCycleMonomers(usePep)
        >>> for pair in mPairs:
        ...     print(pair.monomer1)
        ...     print(pair.monomer2)
        NumberedMonomer(chain=2, number=12, monomer='[E(2,3)]')
        NumberedMonomer(chain=2, number=14, monomer='[E(2,3)]')
        >>> path = BILNParser.GetIntraChainPath(chains[1],
        ...     mPairs[0].monomer1, mPairs[0].monomer2)
        Traceback (most recent call last):
        ...
        ValueError: Cannot match given monomers in given chain.

        Resetting the monomer numbers can fix the problem:
        >>> mPairs = BILNParser.GetIntraChainCycleMonomers(usePep, resetMonomerNumbers=True)
        >>> path = BILNParser.GetIntraChainPath(chains[1],
        ...     mPairs[0].monomer1, mPairs[0].monomer2)
        >>> for monomer in path:
        ...     print(monomer)
        NumberedMonomer(chain=2, number=2, monomer='[E(2,3)]')
        NumberedMonomer(chain=1, number=3, monomer='C')
        NumberedMonomer(chain=2, number=4, monomer='[E(2,3)]')

        Another more practical example:
        >>> usePep = "4xA-C(1,3)-5xA-C(1,3)-P-E-P"
        >>> usePep = BILNParser.DecompressBILN(usePep)
        >>> print(usePep)
        A-A-A-A-C(1,3)-A-A-A-A-A-C(1,3)-P-E-P

        >>> chains = BILNParser.SplitToNumberedChains(usePep)
        >>> mPairs = BILNParser.GetIntraChainCycleMonomers(usePep)
        >>> path = BILNParser.GetIntraChainPath(chains[0],
        ...     mPairs[0].monomer1, mPairs[0].monomer2)
        >>> for monomer in path:
        ...     print(monomer)
        NumberedMonomer(chain=1, number=5, monomer='[C(1,3)]')
        NumberedMonomer(chain=1, number=6, monomer='A')
        NumberedMonomer(chain=1, number=7, monomer='A')
        NumberedMonomer(chain=1, number=8, monomer='A')
        NumberedMonomer(chain=1, number=9, monomer='A')
        NumberedMonomer(chain=1, number=10, monomer='A')
        NumberedMonomer(chain=1, number=11, monomer='[C(1,3)]')

        """
        if not isinstance(numberedChain, NumberedChain):
            raise TypeError(
                "Expecting NumberedChain, got %s" % type(numberedChain))
        # TEMP: skip monomer type check...
        if len(set((numberedChain.chain, monomer1.chain, monomer2.chain))) != 1:
            raise ValueError("Chain IDs do not match.")
        posLow = min(monomer1.number, monomer2.number)
        posHigh = max(monomer1.number, monomer2.number)
        numbered = BILNParser.GetNumberedBILN(numberedChain.BILN,
            returnFormat='monomerobj')
        lowMatch = [m for m in numbered if m.number == posLow]
        highMatch = [m for m in numbered if m.number == posHigh]
        if len(lowMatch) != 1 or len(highMatch) != 1:
            raise ValueError("Cannot match given monomers in given chain.")
        keepResidues = [m for m in numbered if posLow < m.number < posHigh]
        if keepEnds:
            keepResidues += [monomer1, monomer2]
        return tuple(sorted(keepResidues, key=lambda m: m.number))
    ########################################


# end of BILNParser object implementation.
############################################################


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
