import itertools
import numpy as np
import pandas as pd

from scipy import stats
from scipy.spatial import distance
# from sklearn.metrics import *
from sklearn.metrics import mean_squared_error
from sklearn.metrics import calinski_harabasz_score


def _spearman_distance(u, v):
    return 1 - stats.spearmanr(u, v).correlation


def spearman_d(U, v=None, to_numpy=True):
    '''
    Calculate Spearman distance (1 - Spearman correlation coefficient)

    Parameters
    ----------
    U : 1-D or 2-D array or pandas.DataFrame
        If 1-D, you must also provide `v`. If 2-D or pandas.DataFrame, m
        samples (rows) x n columns (features).
    v: None, 1-D array
        If `v` is not None, then both `U` and `v` must be 1-D arrays.
    to_numpy : bool
        If True (default) and `U` is 2-D, return a numpy array. If False,
        return a pandas.DataFrame.

    Returns
    -------
    float, numpy.ndarray, pandas.DataFrame
        If `U` is 2-D, returns `m x m` correlation matrix (pair-wise among
        samples). If `U` is 1-D and `v` is not none, return the Spearman
        distance (float).
    '''
    # The alternate method to calculate Spearman distances is to apply the
    # `_spearman_distance` function to all rows of the array/dataframe like
    # this:
    # >>> dvec = [_spearman_distance(u, v) for u,v in itertools.combinations(U, 2)]
    # >>> dmat = distance.squareform(dvec)
    # This^ is 1-2 orders of magnitude slower than the built-in `.corr` method
    # in pandas. Therefore, if `U` is an array, convert it to a dataframe, and
    # use the `.corr` method. Not sure at what point the performance trade off
    # for converting data type might become a bottleneck.
    if v is None:
        if isinstance(U, np.ndarray):
            U = pd.DataFrame(U)

        if to_numpy:
            return (1 - U.T.corr(method='spearman')).to_numpy()
        else:
            return 1 - U.T.corr(method='spearman')
    else:
        # Assume both `U` and `v` are 1D.
        return _spearman_distance(U, v)


# END
