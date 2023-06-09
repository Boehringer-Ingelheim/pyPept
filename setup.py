"""
Installation script for pyPept.

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
import setuptools

with open("README.md", "r", encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='pyPept',
    version='1.0',
    url='https://github.com/Boehringer-Ingelheim/pyPept',
    license='MIT',
    author='Rodrigo Ochoa, J.B Brown, Thomas Fox',
    author_email='rodrigo.ochoa@boehringer-ingelheim.com, thomas.fox@boehringer-ingelheim.com',
    description='pyPept: a python library to generate atomistic 2D and 3D representations of peptides',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=[
        'numpy >= 1.22.2',
        'pandas >= 1.4.1',
        'requests >= 2.27.1',
        'igraph >= 0.9.10',
        'biopython >= 1.79'
    ],
    include_package_data = True,
    packages=['pyPept','pyPept.data'],
    package_data = {'pyPept': ['data/*']}
)

try:
    import rdkit
except ModuleNotFoundError as exc:
    print("No rdkit installation found! Make sure you install it first!\n" + \
        "Try: 'conda install -c conda-forge rdkit'.")

if rdkit.__version__ < '2019.09.04':
    raise ValueError(
        f"rdkit version must be at least 2019.09.04, but is {rdkit.__version__}.\n" + \
        "Please update first using 'conda update rdkit'")
