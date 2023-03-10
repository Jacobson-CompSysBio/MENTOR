# from functional_partitioning._version import get_version
from functional_partitioning._version import __version__
from functional_partitioning import _cluster as cluster
from functional_partitioning import _metrics as metrics

# __version__ = get_version()
__all__ = [
    '__version__',
    'cluster',
    'metrics'
]

# from . import functional_partitioning
# from . import _datasets as datasets
