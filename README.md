# pyPept

## A python library to generate atomistic 2D and 3D representations of peptides

* From the publication [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2) 
* Journal of Cheminformatics, 2023, 15:79
* Authors: Rodrigo Ochoa, J.B Brown, Thomas Fox

## Purpose

Here we present pyPept, a package to allow the analysis of natural and modified peptides that are assembled based on personalized monomer dictionaries using the Boehringer Ingelheim line notation format (BILN). From the BILN, the peptide construct can then be represented as an RDKit object for further prediction of properties and chemical structures.

## Required third-party tools

The package relies on RDKit (https://rdkit.org/) and the BioPython (https://biopython.org/) python packages to map the BILN peptide and generate the different molecular formats.

## Quick installation

We recommend creating a conda environment with python 3.9
```Bash
conda create -n pypept python=3.9
conda activate pypept
```
The remaining dependencies such as BioPython and Pandas can be installed using the `setup.py` file provided in the code repository, which uses python package managers to easily install the required modules. The script can be called with:
```Bash
   pip install git+ssh://git@github.com:Boehringer-Ingelheim/pyPept.git
```

That's all there is to installation! pyPept can be run using a CLI file provided called `run_pyPept` or by using the modules directly in a python script.
Examples of both cases are described in the next section.

## How to run it

### 1. Using the command-line script provided

The script `run_pyPept` has the following arguments:

```  
usage: run_pyPept.py [-h] (--biln string | --helm string | --fasta string) 
                       [--depiction text] [--prefix text] [--secstruct text] [--noconf] 
                       [--imagesize dim dim] [--logfile filename] [-v]

Generate atomistic 2D and 3D representations of peptides from 
given monomer sequences.

Main arguments:
  -h, --help           show this help message and exit
  --biln string        BILN string with the peptide to analyze.
  --helm string        HELM string with the peptide to analyze.
  --fasta string       FASTA string with the peptide to analyze. 
                       Only natural amino acids are allowed.

Additional options:
  --depiction text     Method to generate the 2D image. 
                       Two options are supported: 'local' (default) or 'rdkit'.
  --prefix text        Name used in the output files. The default is 'peptide'.
  --secstruct text     Use the given secondary structure. 
                       Otherwise, the secondary structure is predicted and used.
  --sdf2D              Generate a 2D SDF file of the peptide.
  --noconf             Do not generate a conformer for the peptide.
  --imagesize dim dim  Image size for 2D depiction, default (1200, 1200).

Logging options:
  --logfile filename   Output messages to given logfile, default is stderr.
  -v, --verbose        Increase output verbosity
```

The only required variable is the peptide, which can be provided directly using the BILN format (--biln), or both HELM (--helm) and FASTA (--fasta) can be provided too. For the latest two, the pipeline script converts the format to a BILN representation. For FASTA only natural amino acids are allowed.

Specifically, pyPept can interconvert between HELM and FASTA formats to BILN, which is internally used by pyPept to facilitate the readability of the molecules and to guarantee a correct chemistry of the peptides during the generation of the RDKit molecular object. Some examples of different modified peptides are shown below:

| BILN | HELM | FASTA |
| --- | --- | --- |
| P-E-P-T-I-D-E | PEPTIDE1{P.E.P.T.I.D.E}$$$$V2.0 |	PEPTIDE |
| ac-D-T-H-F-E-I-A-am | PEPTIDE1{[ac].D.T.H.F.E.I.A.[am]}$$$$V2.0 | None |
| C(1,3)-A-A-A-C(1,3) | PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0 | CAAAC |
| A-G-Q-A-A-K(1,3)-E-F-I-A-A.G-L-E-E(1,3) | PEPTIDE1{A.G.Q.A.A.K.E.F.I.A.A}PEPTIDE2{G.L.E.E}$PEPTIDE1,PEPTIDE2,6:R3-4:R3$$$V2.0 | None |
| N-Iva-F-D-I-meT-N-A-L-W-Y-Aib-K | PEPTIDE1{N.[Iva].F.D.I.[meT].N.A.L.W.Y.[Aib].K}$$$$V2.0 | None |

Additional options can be included based on the description provided in the help menu. An example using the biln sequence 'ac-D-T-H-F-E-I-A-am' is shown (``ac`` and ``am`` represent the terminal acid and amine, respectively, and are defined monomers as part of the pyPept package):

```Bash
run_pyPept --biln 'ac-D-T-H-F-E-I-A-am'
```

### 2. Using the modules individually

If the functions want to be used separately, these are examples for each available class. The first thing is to import the modules and main functions:

**NOTE: In the `examples` folder we provide a set of scripts to run some of the modules for generating 2D and 3D  representations using the peptides shown in the table with different formats (BILN and HELM).**

```Python
# PyPept modules
from pyPept.sequence import Sequence
from pyPept.sequence import correct_pdb_atoms
from pyPept.molecule import Molecule
from pyPept.converter import Converter
from pyPept.conformer import Conformer
from pyPept.conformer import SecStructPredictor

# RDKit modules
from rdkit import Chem
from rdkit.Chem import Draw
```

The biln peptide can be assigned to a `biln` variable in order to create a sequence object:

```Python
# Start the Sequence object
biln = "ac-D-T-H-F-E-I-A-am"
seq = Sequence(biln)
# Correct atom names in the sequence object
seq = correct_pdb_atoms(seq)
```

If the peptide is in HELM notation, it can be converted to BILN using the following function:

```Python
# Call the converter to change from HELM to BILN
from pyPept import Converter

helm = "PEPTIDE1{[ac].D.T.H.F.E.I.A.[am]}$$$$V2.0"
b = Converter(helm=helm)
biln = b.get_biln()
seq = Sequence(biln)
# Correct atom names in the sequence object
seq = correct_pdb_atoms(seq)
```


The Sequence class can receive a `path` variable if the ``data`` folder (including the monomer library) is located in a different location in the system. After creating the sequence, the monomers in the RDKit format can be inspected using a loop as follows:

```Python
# Loop wit the included monomers
mm_list = seq.s_monomers
for i, monomer in enumerate(mm_list):
    mon = monomer['m_romol']
```

The Molecule class can be called after creating the Sequence object. An example to generate the RDKit object, print the SMILES and generate a 2D depiction is as follows:

```Python
# Generate the RDKit object
mol = Molecule(seq)
romol = mol.get_molecule(fmt='ROMol')
print("The SMILES of the peptide is: {}".format(Chem.MolToSmiles(romol)))
Draw.MolToFile(romol, 'peptide.png', size=(1200, 1200))
```

After having the RDKit molecule object, the user can call the Conformer class and generate a PDB file with predicted secondary structure restraints as explained in the paper, using a residue-like format with corrected atom names. The secondary structure can be provided manually if required. The main categories are: B (beta bridge), H (alpha helix). E (beta strand), S (bend),T (turn) and G (3/10 helix)

An example to generate a conformer is shown:

```Python
# Create the peptide conformer with corrected atom names and secondary structure
# Obtain peptide main chain to predict the secondary structure
fasta = Conformer.get_peptide(biln)
secstruct = SecStructPredictor.predict_active_ss(fasta)
# Generate the conformer
romol = Conformer.generate_conformer(romol, secstruct, generate_pdb=True)
```

The RDKit object has now embedded the conformer with the correct atom names. Optionally the user can pass arguments to change the name of the output files, as well as provide a different path to access the data folder.

## Tests

A set of unit tests are available in the `tests` folder. These can be run separately per module by calling each test script, or all can be tested at the same time using the `test.py` file.

```Bash
python test.py
```

## References

If you use pyPept in your work, please cite the following papers:

* [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2) , Journal of Cheminformatics, 2023.
* [BILN – A Human-readable Line Notation for Complex Peptides](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00703), Journal of Chemical Information and Modelling, 2022.

## Support

For inquiries please contact: thomas.fox@boehringer-ingelheim.com .
