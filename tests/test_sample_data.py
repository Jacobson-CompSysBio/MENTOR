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
    'ranks_matrix': 'c59db54a79414e67d733d5e810f0e82a',
    'fullranks': '9af03c1ae884a9f663975bca2498e5cb',
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


# X (features) and y (target).
def test_make_features(random_state=RANDOM_STATE):
    X, y = datasets.make_features()
    assert joblib.hash(X) == CHECKSUMS['features']

# Scores matrix.
def test_make_scores_matrix(random_state=RANDOM_STATE):
    scores_matrix = datasets.make_scores_matrix()
    assert joblib.hash(scores_matrix) == CHECKSUMS['scores_matrix']


# Ranks matrix.
def test_ranks_matrix_na_values():
    ranks_matrix = datasets.make_ranks_matrix()
    count_na_values = ranks_matrix.isna().sum().sum()
    assert count_na_values == 9
    
def test_ranks_matrix_min_rank():
    ranks_matrix = datasets.make_ranks_matrix()
    min_rank = ranks_matrix.min().min()
    assert min_rank == 1

def test_ranks_matrix_max_rank():
    ranks_matrix = datasets.make_ranks_matrix()
    max_rank = ranks_matrix.max().max()
    assert max_rank == 25
    
def test_make_ranks_matrix(random_state=RANDOM_STATE):
    ranks_matrix = datasets.make_ranks_matrix()
    md5sum = joblib.hash(ranks_matrix)
    # print(md5sum)
    assert md5sum == CHECKSUMS['ranks_matrix']


# RWRtoolkit "fullranks" file.
# def test_make_fullranks_table(random_state=RANDOM_STATE):
#     fullranks = datasets.make_fullranks_table()
#     assert joblib.hash(fullranks) == CHECKSUMS['fullranks']
    
class TestFullranksTable(unittest.TestCase):
    def setUp(self):
        # Load test data
        self.data = datasets.make_fullranks_table()

    def test_make_fullranks_table(self):
        self.md5sum = joblib.hash(self.data)
        # print(self.md5sum)
        assert self.md5sum == CHECKSUMS['fullranks']

    def test_score_min(self):
        self.min_score = self.data['Score'].min()
        print(self.min_score)
        testing.assert_almost_equal(0.0, self.min_score, 2)
    
    def test_score_max(self):
        self.max_score = self.data['Score'].max()
        print(self.max_score)
        testing.assert_almost_equal(1.0, self.max_score, 2)
        
    def test_score_na_values(self):
        self.has_na = self.data['Score'].isna().any()
        assert self.has_na != False, 'Fullranks table should not contain NA scores.' 
