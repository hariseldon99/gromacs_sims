"""
protlig_api.steps compatibility wrapper.
"""
try:
    from gmx_protlig.steps import *
except ImportError:
    from ..gmx_protlig.steps import *
