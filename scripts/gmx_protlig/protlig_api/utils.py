"""
protlig_api.utils compatibility wrapper.
"""
try:
    from gmx_protlig.utils import *
except ImportError:
    from ..gmx_protlig.utils import *
