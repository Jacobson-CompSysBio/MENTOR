import itertools
import joblib
import numpy as np
import pandas as pd
import os
import unittest

from numpy import testing
from sklearn import datasets, preprocessing
from matplotlib import pyplot

from mentor import _cluster as cluster
from mentor import _datasets as datasets
from mentor import _metrics as metrics
from mentor import _rwrtoolkit as rwrtoolkit
#from functional_partitioning import _plot as plot

from scipy.spatial import distance
from scipy.cluster import hierarchy


RANDOM_STATE = 42

# Test metrics functions.
def test_spearman_distance_function_equal_pandas_scores():
    '''
    Using the `metrics._spearman_distance` fxn is much slower than
    `pandas.DataFrame.corr`, so I want to check to make sure that the two
    methods are giving (approximately) the same output.
    '''
    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    dvec = [metrics._spearman_distance(u, v) for u,v in itertools.combinations(scores.to_numpy(), 2)]
    scipy_dmat = distance.squareform(dvec)

    pandas_dmat = 1 - scores.T.corr(method='spearman')

    np.testing.assert_array_almost_equal(scipy_dmat, pandas_dmat.values)


def test_spearman_distance_function_equal_pandas_ranks():
    '''
    Using the `metrics._spearman_distance` fxn is much slower than
    `pandas.DataFrame.corr`, so I want to check to make sure that the two
    methods are giving (approximately) the same output.
    '''
    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)

    dvec = [metrics._spearman_distance(u, v) for u,v in itertools.combinations(ranks.to_numpy(), 2)]
    scipy_dmat = distance.squareform(dvec)

    pandas_dmat = 1 - ranks.T.corr(method='spearman')

    np.testing.assert_array_almost_equal(scipy_dmat, pandas_dmat.values)


def test_spearman_distance_matrix_scores():
    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    dmat = metrics.spearman_d(scores)

    assert joblib.hash(dmat) == '181f9b182d74532fdeba6ce50051c324'


def test_spearman_distance_matrix_ranks():
    ranks = datasets.make_scores_matrix()
    ranks = ranks.fillna(0)

    dmat = metrics.spearman_d(ranks)

    assert joblib.hash(dmat) == '17ec19ed3ff1dde6a6f0dfdd210e40ca'


# END.
