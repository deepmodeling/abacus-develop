'''implements the parser for the latest version of ABACUS'''
from typing import List
from ase.atoms import Atoms
def read_abacus_out(fileobj, index=slice(None), results_required=True) \
    -> Atoms | List[Atoms]:
    raise NotImplementedError('The latest version of ABACUS is not supported')