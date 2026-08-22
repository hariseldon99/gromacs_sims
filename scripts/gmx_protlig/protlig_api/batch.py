"""
protlig_api.batch compatibility wrapper.
"""
try:
    from gmx_protlig.batch import *
except ImportError:
    from ..gmx_protlig.batch import *
