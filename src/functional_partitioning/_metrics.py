import itertools
import numpy as np
import pandas as pd
import logging
from scipy import stats
from scipy.spatial import distance
from sklearn.metrics import mean_squared_error
from sklearn.metrics import calinski_harabasz_score
from sklearn.metrics import pairwise_distances

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
    if v is None:
        if isinstance(U, np.ndarray):
            U = pd.DataFrame(U)
        if to_numpy:
            return (1 - U.T.corr(method='spearman')).to_numpy()
        else:
            return 1 - U.T.corr(method='spearman')
    else:
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
    rmse_c = (
        ( (c-1)/(b-1) ) * _root_mean_squared_error(Lc, **kwargs)
    ) + (
        ( (b-c)/(b-1) ) * _root_mean_squared_error(Rc, **kwargs)
    )
    return rmse_c

def get_elbow(Y, min_size=3, **kwargs):
    r'''
    RMSE_{c}={c-1\over b-1}\times RMSE(L_{c})+{b-c\over b-1}\times RMSE(R_{c}) \eqno{\hbox{[1]}}

    Parameters
    ----------
    Y : array-like
        1-D array of values to find the elbow of.
    min_size : int
        Minimum size of the left and right clusters.

    Returns
    -------
    c : int
        The index of the elbow.

    Examples
    --------
    >>> y = get_scores_vs_ranks_curve(scores_matrix)
    >>> max_rank = get_elbow(y)
    '''
    if isinstance(Y, pd.DataFrame):
        raise ValueError('Y must be a numpy array or pandas Series.')
    elif isinstance(Y, pd.Series):
        Y = Y.values
    else:
        pass
    b = len(Y)
    rmse_over_c = []
    for c in range(min_size, b-(min_size+1)):
        rmse_at_c = _root_mean_squared_error_at_c(Y, c, **kwargs)
        LOGGER.debug(f'rmse_at_c: {rmse_at_c}')
        rmse_over_c.append(rmse_at_c)

    idx_of_elbow = int(np.argmin(rmse_over_c) + min_size)
    return idx_of_elbow

def calc_chi(X_true, clusters):
    if isinstance(clusters, (pd.DataFrame, pd.Series)):
        clusters = clusters.to_numpy()
    n_samples = X_true.shape[0]
    n_cuts = clusters.shape[-1]
    chi_scores = []
    for i in range(n_cuts):
        labels_pred = clusters[:, i]
        n_clusters = len(set(labels_pred))
        if n_clusters < 2:
            chi_scores.append(np.nan)
            continue
        elif n_clusters > (n_samples - 1):
            chi_scores.append(np.nan)
            continue
        chi = calinski_harabasz_score(X_true, labels_pred)
        chi_scores.append(chi)
    chi_scores = np.array(chi_scores)
    return chi_scores

def calinski_harabasz_score_(X, labels, **kwargs):
    if labels.ndim == 1:
        pass
    elif labels.ndim == 2:
        if labels.shape[0] == 1 or labels.shape[1] == 1:
            labels = labels.ravel()
    else:
        raise ValueError(f'labels must be 1-D, you gave {labels.shape}')
    return calinski_harabasz_score(X, labels, **kwargs)

def summarize_pairwise_dissimilarities(dissimilarities, labels):
    '''
    Parameters
    ----------
    dissimilarities : array-like
        A condensed distance matrix (1-D array of dissimilarities).
    labels : array-like
        A 1-D array of labels.

    Returns
    -------
    summary : dict
    '''
    if dissimilarities.ndim > 1:
        raise ValueError(f'dissimilarities must be 1-D, you gave {dissimilarities.shape}')
    if labels.ndim == 1:
        pass
    elif labels.ndim == 2:
        if labels.shape[0] == 1 or labels.shape[1] == 1:
            labels = labels.ravel()
    else:
        raise ValueError(f'labels must be 1-D, you gave {labels.shape}')
    within_dissimilarity = []
    between_dissimilarity = []
    for i, pair in enumerate(itertools.combinations(labels, 2)):
        if pair[0] == pair[1]:
            within_dissimilarity.append(dissimilarities[i])
        else:
            between_dissimilarity.append(dissimilarities[i])
    result = {
        'pairwise_dissimilarity_mean': np.mean(dissimilarities),
        'pairwise_dissimilarity_median': np.median(dissimilarities),
        'within_cluster_dissimilarity_mean': np.mean(within_dissimilarity),
        'between_cluster_dissimilarity_mean': np.mean(between_dissimilarity),
        'within_cluster_dissimilarity_median': np.median(within_dissimilarity),
        'between_cluster_dissimilarity_median': np.median(between_dissimilarity)
    }
    return result

def get_scores_vs_ranks_curve(scores, ranks=None):
    '''
    Convert the scores to ranks (or use the provided ranks), and then get the
    mean score at each rank.

    Parameters
    ----------
    scores : array-like
        The scores matrix
    ranks : array-like, optional
        The ranks matrix. If not provided, the scores will be converted to ranks.

    Returns
    -------
    mean_score_vs_rank : array-like
    '''
    stacked_scores = scores.stack()
    if ranks is None:
        stacked_ranks = stacked_scores.groupby(level=0).rank(ascending=False, method='first')
    else:
        stacked_ranks = ranks.stack()
    mean_score_vs_rank = stacked_scores.groupby(stacked_ranks).mean()
    return mean_score_vs_rank
