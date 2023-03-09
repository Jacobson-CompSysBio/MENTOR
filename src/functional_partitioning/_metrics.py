import itertools
import numpy as np
import pandas as pd
import logging

from scipy import stats
from scipy.spatial import distance
# from sklearn.metrics import *
from sklearn.metrics import mean_squared_error
from sklearn.metrics import calinski_harabasz_score

LOGGER = logging.getLogger(__name__)

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


def _root_mean_squared_error(y_true, y_pred=None, **kwargs):
    '''
    y_true, y_pred, *, sample_weight=None, multioutput='uniform_average', squared=True
    '''
    if y_pred is None:
        y_pred = np.linspace(y_true[0], y_true[-1], len(y_true))
    if len(y_true) != len(y_pred):
        raise ValueError('y_true and y_pred must be the same length.')
    mse = mean_squared_error(y_true, y_pred, **kwargs)
    rmse = np.sqrt(mse)
    return rmse


def _root_mean_squared_error_at_c(y_true, c, **kwargs):
    '''
    Root mean squared error at c
    '''
    b = len(y_true)
    l_start = 0
    l_stop = c+1
    r_start = c
    r_stop = b+1
    Lc = y_true[l_start:l_stop]
    Rc = y_true[r_start:r_stop]
    LOGGER.debug(rf'c={c}; Lc=[{l_start}:{l_stop}] ({len(Lc)}); Rc=[{r_start}:{r_stop}] ({len(Rc)})')
    # RMSE at c is the sum of RMSE to the left and RMSE to the right.
    rmse_c = (
        ( (c-1)/(b-1) ) * _root_mean_squared_error(Lc, **kwargs)
    ) + (
        ( (b-c)/(b-1) ) * _root_mean_squared_error(Rc, **kwargs)
    )
    return rmse_c


def get_elbow(y_true, min_size=3, **kwargs):
    r'''
    RMSE_{c}={c-1\over b-1}\times RMSE(L_{c})+{b-c\over b-1}\times RMSE(R_{c}) \eqno{\hbox{[1]}}

    Parameters
    ----------
    y_true : array-like
        The true values.
    min_size : int
        Minimum size of the left and right clusters.

    Returns
    -------
    c : int
        The index of the elbow.

    Examples
    --------
    >>> y_true = scores_matrix.values
    >>> c = get_elbow(y_true)
    '''

    if isinstance(y_true, pd.DataFrame):
        raise ValueError('y_true must be a numpy array or pandas Series.')
    elif isinstance(y_true, pd.Series):
        y_true = y_true.values
    else:
        pass

    b = len(y_true)

    rmse_over_c = []

    for c in range(min_size, b-(min_size+1)):
        rmse_at_c = _root_mean_squared_error_at_c(y_true, c, **kwargs)
        rmse_over_c.append(rmse_at_c)
    # Adjust index by min_size.
    idx_of_elbow = int(np.argmin(rmse_over_c) + min_size)
    return idx_of_elbow


def calc_chi(X_true, clusters):
    if isinstance(clusters, (pd.DataFrame, pd.Series)):
        clusters = clusters.to_numpy()
    # CHI is only valid for clusterings with n-clusters between 2 and n samples-1.
    # Filling with NaN is more accurate, but plotting the values is misleading if the user
    # is unaware that the missing values are not plotted.
    n_samples = X_true.shape[0]
    n_cuts = clusters.shape[-1]
    chi_scores = []
    for i in range(n_cuts):
        labels_pred = clusters[:, i]
        n_clusters = len(set(labels_pred))
        if n_clusters < 2:
            chi_scores.append(np.nan)
            # chi_scores.append(0)
            continue
        elif n_clusters > (n_samples - 1):
            chi_scores.append(np.nan)
            # chi_scores.append(0)
            continue
        chi = calinski_harabasz_score(X_true, labels_pred)

        chi_scores.append(chi)
    chi_scores = np.array(chi_scores)
    return chi_scores

# END
