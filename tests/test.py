"""
Script to run all package tests in one place.

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

import sys
from unittest import TestLoader, TestResult
from pathlib import Path

##########################################################################
# Function
##########################################################################

def run_tests():
    """
    Function to run all the available unittests
    """

    test_loader = TestLoader()
    test_result = TestResult()

    test_directory = str(".")

    test_suite = test_loader.discover(test_directory, pattern='*_Test.py')
    test_suite.run(result=test_result)

    if test_result.wasSuccessful():
        print(f"\nAll {test_result.testsRun} groups of tests " +
            "(conformers, molecules, seq->mol, sequence perception) " +
            "were successful!")
        sys.exit(0)
    else:
        print(f"\n{len(test_result.failures)} failures:")
        for i,fail in enumerate(test_result.failures):
            print(f"{i+1}. {fail[0]}\n{fail[1]}")
        if test_result.errors:
            print(f"\nErrors: {test_result.errors}")
        sys.exit(-1)

if __name__ == '__main__':
    run_tests()
