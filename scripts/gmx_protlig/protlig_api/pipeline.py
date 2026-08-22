"""
protlig_api.pipeline compatibility wrapper.
"""
try:
    from gmx_protlig.pipeline import *
except ImportError:
    from ..gmx_protlig.pipeline import *
