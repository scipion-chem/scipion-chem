# Lightweight shim so scripts can import `constants.descriptor_categories`
# without importing the whole `pwchem` package (which needs pyworkflow).
from importlib.machinery import SourceFileLoader
import os

_repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
_src = os.path.join(_repo_root, 'constants.py')
_pwchem_constants = SourceFileLoader('pwchem_constants', _src).load_module()

# Expose the lowercase name expected by ligand_descriptor_calc.py
descriptor_categories = getattr(_pwchem_constants, 'DESCRIPTOR_CATEGORIES', {})
