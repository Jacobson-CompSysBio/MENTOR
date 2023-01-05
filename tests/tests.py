import joblib
import itertools
import numpy as np

from sklearn import datasets, preprocessing
from matplotlib import pyplot

from functional_partitioning import datasets, cluster, metrics, rwrtoolkit
from functional_partitioning import functional_partitioning as fp

from scipy.spatial import distance
from scipy.cluster import hierarchy


CHECKSUMS = {
    'features': '05f5c7599ac919bddb182aec8d2aaf16',
    'scores_matrix': '0bc8e7374e6aa2b56ca100eabad48e4b',
    'ranks_matrix': '946e98267a27f4d837bdf0cbb4d1a529',
    'fullranks': 'e8646ece5873b493e38a71f877372e66',
    'linkage_matrix_scores_euclidean_average': 'b1dd3e10cc1029f1a57202f4be3c9cfa',
    'linkage_matrix_ranks_euclidean_average': '477efbcf8ed2e66e4e2ea29d77d57abc',
    'spearman_distance_matrix_scores': '181f9b182d74532fdeba6ce50051c324',
    'spearman_distance_matrix_ranks': '17ec19ed3ff1dde6a6f0dfdd210e40ca',
    # Clusters based on Spearman distances from the `scores` and `ranks` matrices should be the same.
    'spearman_distance_clusters': 'c35db22b2b15458dfe3817f7781dd4ce',
    # RWRtoolkit wrapper.
    'fullranks_to_matrix_X_ranks': '9101f73917fc9678808ce47d720ec411',
    'fullranks_to_matrix_X_scores' : '3c6d835ec34107f53bd13f9f1e1a16ce',
    'fullranks_to_matrix_labels' : 'bdcc8bd8b0edc2e52b82d39b282493b7',
    'fullranks_to_matrix_max_rank' : 'c763ff8f88aacdbe4b423f930ff3afeb'
}


RANDOM_STATE = 42


# Test data sets.


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


# Test linkage matrices.


def test_scipy_linkage_matrix_scores():
    metric = 'euclidean'
    method = 'average'

    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    dvec = distance.pdist(scores, metric=metric)
    linkage_matrix = hierarchy.linkage(dvec, method=method)

    assert joblib.hash(linkage_matrix) == CHECKSUMS['linkage_matrix_scores_euclidean_average']


def test_scipy_linkage_matrix_ranks():
    metric = 'euclidean'
    method = 'average'

    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)

    dvec = distance.pdist(ranks, metric=metric)
    linkage_matrix = hierarchy.linkage(dvec, method=method)

    assert joblib.hash(linkage_matrix) == CHECKSUMS['linkage_matrix_ranks_euclidean_average']


def test_hierarchicalclustering_linkage_matrix_scores():
    metric = 'euclidean'
    method = 'average'

    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    mod = cluster.HierarchicalClustering(affinity='euclidean')
    mod.fit(scores)

    assert joblib.hash(mod.linkage_matrix) == CHECKSUMS['linkage_matrix_scores_euclidean_average']


def test_hierarchicalclustering_linkage_matrix_ranks():
    metric = 'euclidean'
    method = 'average'

    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)

    mod = cluster.HierarchicalClustering(affinity='euclidean')
    mod.fit(ranks)

    assert joblib.hash(mod.linkage_matrix) == CHECKSUMS['linkage_matrix_ranks_euclidean_average']


def test_hierarchicalclustering_clusters_scores():
    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    mod = cluster.HierarchicalClustering(affinity=metrics.spearman_d)
    mod.fit(scores)

    clusters = hierarchy.cut_tree(mod.linkage_matrix)

    assert joblib.hash(clusters) == CHECKSUMS['spearman_distance_clusters']


def test_hierarchicalclustering_clusters_ranks():
    ranks = datasets.make_scores_matrix()
    ranks = ranks.fillna(ranks.max(axis=1))

    mod = cluster.HierarchicalClustering(affinity=metrics.spearman_d)
    mod.fit(ranks)

    clusters = hierarchy.cut_tree(mod.linkage_matrix)

    assert joblib.hash(clusters) == CHECKSUMS['spearman_distance_clusters']


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

    assert joblib.hash(dmat) == CHECKSUMS['spearman_distance_matrix_scores']


def test_spearman_distance_matrix_ranks():
    ranks = datasets.make_scores_matrix()
    ranks = ranks.fillna(0)

    dmat = metrics.spearman_d(ranks)

    assert joblib.hash(dmat) == CHECKSUMS['spearman_distance_matrix_ranks']


# Test RWRtoolkit wrapper functions. Minimal testing here; just make sure the
# fullranks table is properly converted.


def test_rwrtoolkit_fullranks_to_matrix_X_ranks():
    fullranks = datasets.make_fullranks_table()
    X_ranks, X_scores, labels, max_rank = rwrtoolkit.fullranks_to_matrix(fullranks)
    assert joblib.hash(X_ranks) == CHECKSUMS['fullranks_to_matrix_X_ranks']


def test_rwrtoolkit_fullranks_to_matrix_X_scores():
    fullranks = datasets.make_fullranks_table()
    X_ranks, X_scores, labels, max_rank = rwrtoolkit.fullranks_to_matrix(fullranks)
    assert joblib.hash(X_scores) == CHECKSUMS['fullranks_to_matrix_X_scores']


def test_rwrtoolkit_fullranks_to_matrix_labels():
    fullranks = datasets.make_fullranks_table()
    X_ranks, X_scores, labels, max_rank = rwrtoolkit.fullranks_to_matrix(fullranks)
    assert joblib.hash(labels) == CHECKSUMS['fullranks_to_matrix_labels']


def test_rwrtoolkit_fullranks_to_matrix_max_rank():
    fullranks = datasets.make_fullranks_table()
    X_ranks, X_scores, labels, max_rank = rwrtoolkit.fullranks_to_matrix(fullranks)
    assert joblib.hash(max_rank) == CHECKSUMS['fullranks_to_matrix_max_rank']


# END.
