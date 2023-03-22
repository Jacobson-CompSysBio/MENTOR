import itertools
import joblib
import numpy as np
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

CHECKSUMS = {
    'features': '05f5c7599ac919bddb182aec8d2aaf16',
    'scores_matrix': '0bc8e7374e6aa2b56ca100eabad48e4b',
    'ranks_matrix': '946e98267a27f4d837bdf0cbb4d1a529',
    'fullranks': 'e8646ece5873b493e38a71f877372e66',
    'linkage_matrix_scores_euclidean_average': '448a5cc77be84dddfe28786308afb4b9',
    'linkage_matrix_ranks_euclidean_average': '477efbcf8ed2e66e4e2ea29d77d57abc',
    'linkage_matrix_scores_spearman_average' : 'cc82aa17649c4ff0e7faa8c4288b8185',
    'linkage_matrix_ranks_spearman_average' : 'b7c48b5580dae5f4fc066d229f098873',
    'spearman_distance_matrix_scores': '181f9b182d74532fdeba6ce50051c324',
    'spearman_distance_matrix_ranks': '17ec19ed3ff1dde6a6f0dfdd210e40ca',
    # Clusters based on Spearman distances from the `scores` and `ranks` matrices should be the same.
    'spearman_distance_clusters': 'c35db22b2b15458dfe3817f7781dd4ce',
    # RWRtoolkit wrapper.
    'fullranks_to_matrix_X_ranks': '9101f73917fc9678808ce47d720ec411',
    'fullranks_to_matrix_X_scores' : '3c6d835ec34107f53bd13f9f1e1a16ce',
    'fullranks_to_matrix_labels' : 'bdcc8bd8b0edc2e52b82d39b282493b7',
    # Helper functions.
    'calc_threshold_mean' : '6b88a48bc7f8bcf7ec27d6e8592fae6a',
    'calc_threshold_best_chi' : 'fb42bd8c5d5c3bec691d6c7c1d37bede',
    'get_clusters_labels_none_threshold_05' : 'a86b536ff791b21ab57fe612a93021df',
    'get_clusters_labels_default_threshold_05' : 'cee8925663a74cb41bd759e6291b8e6d',
    'label_mapper_labels_default' : '38fe01a29198ab24055b1faf53d68cc5'
}



# Test linkage matrices. The linkage matrices returned from
# `scipy.hierarchy.linkage` and `HierarchicalClustering` should be the same.

def test_scipy_linkage_matrix_scores():
    metric = 'euclidean'
    method = 'average'
    decimals = 13

    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    dvec = distance.pdist(scores, metric=metric)
    linkage_matrix = hierarchy.linkage(dvec, method=method)

    # The linkage matrix directly from scipy and the one returned from
    # HierarchicalClustering are the same to 13 decimal places.
    linkage_matrix = linkage_matrix.round(decimals=decimals)
    md5sum = joblib.hash(linkage_matrix)

    assert md5sum == CHECKSUMS['linkage_matrix_scores_euclidean_average']


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
    decimals = 13

    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    mod = cluster.HierarchicalClustering(metric='euclidean')
    mod.fit(scores)

    # The linkage matrix directly from scipy and the one returned from
    # HierarchicalClustering are the same to 13 decimal places.
    linkage_matrix = mod.linkage_matrix.round(decimals=decimals)
    md5sum = joblib.hash(linkage_matrix)

    assert md5sum == CHECKSUMS['linkage_matrix_scores_euclidean_average']


def test_hierarchicalclustering_linkage_matrix_ranks():
    metric = 'euclidean'
    method = 'average'

    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)

    mod = cluster.HierarchicalClustering(metric='euclidean')
    mod.fit(ranks)

    assert joblib.hash(mod.linkage_matrix) == CHECKSUMS['linkage_matrix_ranks_euclidean_average']


def test_hierarchicalclustering_clusters_scores():
    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    mod = cluster.HierarchicalClustering(metric=metrics.spearman_d)
    mod.fit(scores)

    clusters = hierarchy.cut_tree(mod.linkage_matrix)

    assert joblib.hash(clusters) == CHECKSUMS['spearman_distance_clusters']


def test_hierarchicalclustering_clusters_ranks():
    ranks = datasets.make_scores_matrix()
    ranks = ranks.fillna(ranks.max(axis=1))

    mod = cluster.HierarchicalClustering(metric=metrics.spearman_d)
    mod.fit(ranks)

    clusters = hierarchy.cut_tree(mod.linkage_matrix)

    assert joblib.hash(clusters) == CHECKSUMS['spearman_distance_clusters']



def test_cluster_hierarchical_scores_spearman_average():
    scores = datasets.make_scores_matrix()
    scores = scores.fillna(scores.max(axis=1))

    linkage_matrix = cluster.cluster_hierarchical(
        scores,
        corr_method='spearman',
        linkage_method='average'
    )

    assert joblib.hash(linkage_matrix) == CHECKSUMS['linkage_matrix_scores_spearman_average']


def test_cluster_hierarchical_ranks_spearman_average():
    ranks = datasets.make_ranks_matrix()
    ranks = ranks.fillna(0)

    linkage_matrix = cluster.cluster_hierarchical(
        ranks,
        corr_method='spearman',
        linkage_method='average'
    )

    assert joblib.hash(linkage_matrix) == CHECKSUMS['linkage_matrix_ranks_spearman_average']


