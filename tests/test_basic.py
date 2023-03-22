import itertools
import joblib
import numpy as np
import pandas as pd
import os
import unittest

from numpy import testing
from sklearn import datasets, preprocessing
from matplotlib import pyplot

from functional_partitioning import _cluster as cluster
from functional_partitioning import _datasets as datasets
from functional_partitioning import _metrics as metrics
from functional_partitioning import _rwrtoolkit as rwrtoolkit
from functional_partitioning import _plot as plot

from scipy.spatial import distance
from scipy.cluster import hierarchy


RANDOM_STATE = 42

# Test entry point
def test_entrypoint():
    exit_status = os.system('functional_partitioning --help')
    assert exit_status == 0

def test_entrypoint():
    exit_status = os.system('functional_partitioning --version')
    assert exit_status == 0


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


# Test helper functions.

def test_make_label_mapper_labels_default():
    nodetable = pd.DataFrame(
        data=[2, 2, 1, 2, 1, 1, 1, 2, 1],
        index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'],
        columns=['cluster']
    )
    label_mapper = plot.make_label_mapper(
        nodetable=nodetable,
        use_locs=[0, 1], # List.
    )
    md5sum = joblib.hash(label_mapper)
    assert  md5sum == '77cbc5e1b7a83faa49a017350c891987'


# END.
