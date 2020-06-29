"""
I Shrunk me Molecule!
A module to shrink a molecule/ligand during MD simulations in OpenMM
"""

# Add imports here
from .mol_shrink_ray import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
