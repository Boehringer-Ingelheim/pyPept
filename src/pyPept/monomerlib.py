"""
A class to parse and manipulate BILN sequences.

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023

Updated 2025
From presentation: XXX
Certara User Group Meeting, Heidelberg, 2025
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


########################################################################################
# Modules
########################################################################################

# System libraries
import logging
from importlib.resources import files

# Third-party libraries
from rdkit.Chem import PandasTools


##########################################################################
# Functions and classes
##########################################################################

############################################################
class MonomerConstants:
    """
    A container class to hold defaults values related to loading of monomers 
    for a monomer library.
    """
    def_path = "pyPept.data"
    def_lib_filename = "monomers.sdf"
    _sdf_change_separator = ","

    attr_monomer_symbol = "m_abbr"
    attr_linkable_R_groups = "m_Rgroups"
    attr_linkable_R_grop_idx = "m_RgroupIdx"
    attr_linkable_R_group_attachidx = "m_attachmentPointIdx"

############################################################

############################################################
class MonomerLibrary:
    """
    A class that contains functionality related to loading or forming a monomer
    library and querying it for information.
    """

    ############################################################
    def __init__(self, 
                data_dir=MonomerConstants.def_path,
                monomer_lib_file=MonomerConstants.def_lib_filename,
                logger=logging.Logger("dummy")):
        """
        Instantiate a MonomerLibrary object by loading from a file.

        :param data_dir: the directory location of a monomer library.
        :type data_dir: str
        :param monomer_lib_file: the SDF with monomers within the 
            directory to search through.
        :type monomer_lib_file: str
        """

        monomer_df_filepath = files(data_dir).joinpath(monomer_lib_file)

        if monomer_df_filepath.is_file() is False:
            backup_path = files(MonomerConstants.def_path).joinpath(
                                MonomerConstants.def_lib_filename)
            logger.warning(
                "Loading monomer library, cannot resolve %s, trying %s ." % (
                    monomer_df_filepath, backup_path))
            if backup_path.is_file() is False:
                logger.warning(
                    "Loading monomer library, cannot %s . Stopping." % (
                        backup_path))
                raise RuntimeError("Failed to find monomer library.")
            else: # Backup was still available
                monomer_df_filepath = backup_path

        self.monomer_df = self._load_monomer_sdf(str(monomer_df_filepath))
    ############################################################
    
    ############################################################
    def _load_monomer_sdf(self, path):
        """
        Read a monomer SDF file and store it as a Pandas dataframe object.

        :param path: OS-dependent path of an SDF file containing monomers.
        :type path: str

        :return: monomer dictionary as a dataframe
        """
        df_group = PandasTools.LoadSDF(path)

        groups = [
            MonomerConstants.attr_linkable_R_groups,
            MonomerConstants.attr_linkable_R_grop_idx,
            MonomerConstants.attr_linkable_R_group_attachidx,
            ]

        for idx in df_group.index:
            for group in groups:
                change = df_group[group][idx].split(
                    MonomerConstants._sdf_change_separator)
                if group == MonomerConstants.attr_linkable_R_groups:
                    updated_change = [
                        None if v == 'None' else v for v in change]
                else:
                    updated_change = [
                        None if v == 'None' else int(v) for v in change]
                df_group.loc[idx, group] = updated_change
        df_group = df_group.set_index('symbol')
        df_group = df_group.rename(columns={"ROMol": "m_romol"})

        return df_group
    ############################################################

    ############################################################
    def __getitem__(self, library_attribute):
        """
        Retrieve the monomer dictionary attribute requested.
        Currently the monomer dictionary uses a DataFrame as its underlying
        implementation, and this method is equivalent to getting a column of
        values.

        :param library_attribute: the attribute of the monomer dictionary to
                                  retrieve/fetch.
        :type library_attribute: str
        """
        return self.monomer_df[library_attribute]
    ############################################################
    
    ############################################################
    def GetMonomerNames(self):
        """Return the list of monomer names in this library.
        
        :return: a tuple of str
        """
        return tuple(self[MonomerConstants.attr_monomer_symbol])
    ############################################################

    ############################################################
    def GetRGroups(self, monomer):
        """Return the valid R-group identifiers for a given monomer.
        
        :param monomer: a monomer symbol, such as A, DAla, etc.
        :type monomer: str
        :return: a tuple of strings "R1", "R2", ...

        >>> lib = MonomerLibrary()
        >>> lib.GetRGroups("A")
        ('R1', 'R2')
        >>> lib.GetRGroups("D")
        ('R1', 'R2', 'R3')
        >>> lib.GetRGroups("K")
        ('R1', 'R2', 'R3')
        """
        if monomer not in self.GetMonomerNames():
            raise ValueError(
                "%s not valid in this monomer dictionary." % monomer)
        validR = list()
        for idx, group in enumerate(
                self.monomer_df.loc[monomer][
                    MonomerConstants.attr_linkable_R_groups], start=1):
            if group is not None:
                validR.append("R%i" % idx)
        return tuple(validR)
    ############################################################

## end of MonomerLibrary class ##############################

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

############################################################
# End of monomerlib.py
############################################################
