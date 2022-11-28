#!/usr/bin/env python
# coding: utf-8

import joblib

from sklearn import datasets, preprocessing
from matplotlib import pyplot

from functional_partitioning import functional_partitioning as fp
from functional_partitioning import datasets, cluster


CHECKSUMS = {
    'features': '05f5c7599ac919bddb182aec8d2aaf16',
    'scores_matrix': '0bc8e7374e6aa2b56ca100eabad48e4b',
    'ranks_matrix': '946e98267a27f4d837bdf0cbb4d1a529',
    'fullranks': 'e8646ece5873b493e38a71f877372e66',
    'linkage_matrix_scores_euclidean_average': 'b1dd3e10cc1029f1a57202f4be3c9cfa',
    'linkage_matrix_ranks_euclidean_average': '477efbcf8ed2e66e4e2ea29d77d57abc'
}


RANDOM_STATE = 42

# Data sets.
def test_make_features(random_state=RANDOM_STATE):
    X, y = datasets.make_features()
    assert joblib.hash(X) == CHECKSUMS['features']


def test_make_scores_matrix(random_state=RANDOM_STATE):
    scores_matrix = datasets.make_scores_matrix()
    assert joblib.hash(scores_matrix) == CHECKSUMS['scores_matrix']


def test_make_ranks_matrix(random_state=RANDOM_STATE):
    ranks_matrix = datasets.make_ranks_matrix()
    assert joblib.hash(ranks_matrix) == CHECKSUMS['ranks_matrix']


def test_make_fullranks_table(random_state=RANDOM_STATE):
    fullranks = datasets.make_fullranks_table()
    assert joblib.hash(fullranks) == CHECKSUMS['fullranks']


# Linkage matrices.
def test_scipy_linkage_matrix_scores():
    from scipy.spatial import distance
    from scipy.cluster import hierarchy

    metric = 'euclidean'
    method = 'average'

    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    dvec = distance.pdist(scores, metric=metric)
    linkage_matrix = hierarchy.linkage(dvec, method=method)

    assert joblib.hash(linkage_matrix) == CHECKSUMS['linkage_matrix_scores_euclidean_average']


def test_scipy_linkage_matrix_ranks():
    from scipy.spatial import distance
    from scipy.cluster import hierarchy

    metric = 'euclidean'
    method = 'average'

    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)

    dvec = distance.pdist(ranks, metric=metric)
    linkage_matrix = hierarchy.linkage(dvec, method=method)

    assert joblib.hash(linkage_matrix) == CHECKSUMS['linkage_matrix_ranks_euclidean_average']


def test_hierarchicalclustering_linkage_matrix_scores():
    from scipy.spatial import distance
    from scipy.cluster import hierarchy

    metric = 'euclidean'
    method = 'average'

    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))
    
    mod = cluster.HierarchicalClustering(affinity='euclidean')
    mod.fit(scores)

    assert joblib.hash(mod.linkage_matrix) == CHECKSUMS['linkage_matrix_scores_euclidean_average']


def test_hierarchicalclustering_linkage_matrix_ranks():
    from scipy.spatial import distance
    from scipy.cluster import hierarchy

    metric = 'euclidean'
    method = 'average'

    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)
    
    mod = cluster.HierarchicalClustering(affinity='euclidean')
    mod.fit(ranks)

    assert joblib.hash(mod.linkage_matrix) == CHECKSUMS['linkage_matrix_ranks_euclidean_average']


# END.
