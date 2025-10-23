"""
A place where basic abstract types related to BILNs are defined.
"""

# System libraries needed by this module.
from collections import namedtuple

# Third-party libraries needed by this module, e.g. numpy.

# Project-specific modules additionally needed.


# ----- Begin code for this module. -----

# Named tuple classes
NumberedMonomer = namedtuple('NumberedMonomer', ('chain', 'number', 'monomer'))
TitledBILNExample = namedtuple('TitledBILNExample', ('title', 'BILN'))
NumberedChain = namedtuple('NumberedChain', ('chain', 'BILN'))
HLELinker = namedtuple('HLELinker', ('chain', 'BILN'))
HLE = namedtuple('HLE', ('chain', 'branchmonomer', 'linkers', 'FA', 'BILN'))
MonomerPair = namedtuple('MonomerPair', ('monomer1', 'monomer2'))

# Error classes that need to be raised.
class BILNSequenceError(ValueError): pass
class InvalidMonomerName(BILNSequenceError): pass
class AmbiguousBranchIDs(BILNSequenceError): pass

class InvalidMonomerRGroup(BILNSequenceError): pass
class InvalidChainBeginRGroup(InvalidMonomerRGroup): pass
class InvalidChainEndRGroup(InvalidMonomerRGroup): pass
class InvalidChainMiddleRGroup(InvalidMonomerRGroup): pass
class RepeatRGroupMonomer(BILNSequenceError): pass

class BILNMultiError(BILNSequenceError): pass


## Verify that the module's interfaces work as the doctests demonstrate.
if __name__ == "__main__":
  
    import doctest, os, sys, __main__
    numFail, numTests = doctest.testmod()
    modulePath = os.path.abspath(__main__.__file__)
    if numFail > 0:
        print('%s : Expected functionality in doctests fail!' % modulePath)
        sys.exit(-1)
    elif numTests > 0:
        print("%s : All %i doctests passed." % (modulePath, numTests))
        sys.exit(0)
    else:
        print("%s : WARNING - No doctests defined." % modulePath)
        sys.exit(1)
