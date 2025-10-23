#!/usr/bin/env python3
"""
Check the validity of BILN strings, using different criteria.
"""

# System libraries needed by this module.
import argparse  # If used as standalone application.
import logging   # For messages and debugging.
import __main__
import sys

# Third-party libraries needed to execute CLI, e.g. numpy.

# Project-specific modules additionally needed.
from pyPept.biln import BILNConstants, BILNParser

# ----- Begin code for this module. -----
_defColMolID = 1
_defColBILN = 2

############################################################
def getRequiredInputsParser():
    """Constructs parser for absolute minimum required inputs."""

    parser = argparse.ArgumentParser(add_help=False)
    requiredArgs = parser.add_argument_group('Required arguments')

    inputTypeGroup = parser.add_mutually_exclusive_group(required=True)
    inputTypeGroup.add_argument(
        '--biln', type=str, metavar='text',
        required=False,
        help="BILN string to check for validity.")
    inputTypeGroup.add_argument(
        '--table', type=str, metavar='filename',
        required=False,
        help="Table containing BILN strings and molecule IDs.")
    return parser
############################################################


############################################################
def getOptionalInputsParser():
    """Constructs parser for optional arguments."""

    parser = argparse.ArgumentParser(add_help=False)
    options = parser.add_argument_group('Optional arguments')
    options.add_argument(
        '--title', type=str,
        required=False, default="NONE",
        help="For single BILN input, molecule title to assign. Default NONE.")
    
    tableCols = parser.add_argument_group('Table-type input')
    tableCols.add_argument(
        '-d', '--delim', type=str, metavar='char',
        required=False, default="\t",
        help="Delimiter of table, default tab.")
    tableCols.add_argument(
        '--colmolid', type=int,
        required=False, default=_defColMolID,
        help="Column with molecule ID, default %i." % _defColMolID)
    tableCols.add_argument(
        '--colbiln', type=int,
        required=False, default=_defColBILN,
        help="Column with BILN, default %i." % _defColBILN)
    tableCols.add_argument(
        '--noheader', action='store_false', dest='header',
        help="Table has no header, and begin BILN conversion from first line.")
    
    # A repeated-use log option parser.
    logOptions = parser.add_argument_group('Logging options')
    logOptions.add_argument(
        '--logfile', type=str, metavar='filename',
        required=False, default=None,
        help="Output messages to given logfile, default is stderr.")
    logOptions.add_argument(
        "-v", "--verbose", action="store_true", 
        help="Increase output verbosity")

    return parser
############################################################


############################################################
def getStandaloneParser():
    """Constructs parser to run this script as standalone application."""

    parser = argparse.ArgumentParser(
        description="""Check BILN strings for their validity.""",
        parents=(
            getRequiredInputsParser(), getOptionalInputsParser(),
            ))

    return parser
############################################################

############################################################
def commonPrint(title: str, validity: bool, valid_branching: bool,
        invalid_monomers: tuple, valid_monomers: tuple,
        out_delimiter="\t"):
    """
    A print function that uniformly handles output of this tool.
    """
    print(out_delimiter.join([
        title, str(validity), str(valid_branching),
        " ".join(sorted(set(invalid_monomers))) or "None", 
        " ".join(sorted(set(valid_monomers))) or "None",
        ]))
############################################################


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runAsStandalone():
    """Contains script execution flow as main standalone application.
    """

    useParser = getStandaloneParser()
    args = useParser.parse_args()

    # Setup typical logger for messages to stderr.
    logStream = logging.StreamHandler(
        open(args.logfile, "w") if args.logfile else sys.stderr)
    logStream.setFormatter(logging.Formatter(
        "%(asctime)-10s %(levelname)s:%(message)s",
        datefmt="%H:%M:%S"))
    logger = logging.Logger(name="BILN_validation")
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    logger.addHandler(logStream)
    logger.debug("Invocation arguments: %s" % args)

    ########################################
    ## Begin main driver execution #########
    ########################################
    print("\t".join(("MolTitle", "Valid?",
        "ValidBranchAnnotations?", "InvalidMonomers", "ValidMonomers")))

    # Workflow
    if args.biln:
        valid = BILNParser.IsValidSequence(args.biln)
        valMonomers = tuple([m.monomer for m in BILNParser.GetNumberedBILN(
                                        args.biln, returnFormat='monomerobj')
                        if BILNParser.IsValidMonomer(m.monomer)])
        valBranch = BILNParser.HasValidBranchAnnotations(args.biln)
        invalMonomers = BILNParser.GetInvalidMonomers(args.biln)
        commonPrint(args.title, valid, valBranch, invalMonomers, valMonomers)

    # Operating mode 2 : BILNs in a tabular input
    elif args.table:
        with open(args.table) as inTable:
            if args.header:
                inTable.readline()
            for line in inTable:
                try:
                    tokens = line.strip().split(args.delim)
                    title = tokens[args.colmolid - 1]
                    biln = tokens[args.colbiln - 1]
                    valid = BILNParser.IsValidSequence(biln)
                    valMonomers = tuple([m.monomer for m in 
                                    BILNParser.GetNumberedBILN(
                                        biln, returnFormat='monomerobj')
                                    if BILNParser.IsValidMonomer(m.monomer)])
                    valBranch = BILNParser.HasValidBranchAnnotations(biln)
                    invalMonomers = BILNParser.GetInvalidMonomers(biln)
                    commonPrint(title, valid, valBranch,
                        invalMonomers, valMonomers)
                except Exception as e:
                    logger.warning("Failed to process %s: %s" % (tokens, e))

    # Termination and cleanup.
    logger.info("Successful completion of %s." % __main__.__file__)
    
# end of running script as standalone application.


# What to do if script is run as a main driver program.
if __name__ == "__main__":

    runAsStandalone()

# ----- End code for this module. -----
