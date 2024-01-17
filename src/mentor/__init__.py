from mentor._version import __version__
from mentor import _cluster as cluster
from mentor import _metrics as metrics

# __version__ = get_version()
__all__ = [
    '__version__',
    'cluster',
    'metrics'
]
