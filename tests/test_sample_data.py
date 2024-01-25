import itertools
import joblib
import numpy as np
import os
import unittest

from numpy import testing
from sklearn import datasets, preprocessing
#from matplotlib import pyplot

from mentor import _cluster as cluster
from mentor import _datasets as datasets
from mentor import _metrics as metrics
from mentor import _rwrtoolkit as rwrtoolkit
from mentor import _plot as plot

from scipy.spatial import distance
from scipy.cluster import hierarchy

RANDOM_STATE = 42

# X (features) and y (target).
def test_make_features(random_state=RANDOM_STATE):
    X, y = datasets.make_features()
    assert joblib.hash(X) == '05f5c7599ac919bddb182aec8d2aaf16'

# Scores matrix.
def test_make_scores_matrix(random_state=RANDOM_STATE):
    scores_matrix = datasets.make_scores_matrix()
    assert joblib.hash(scores_matrix) == '0bc8e7374e6aa2b56ca100eabad48e4b'


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
    assert md5sum == 'c59db54a79414e67d733d5e810f0e82a'


class TestMakeFullranksTable(unittest.TestCase):
    def setUp(self):
        # Load test data
        self.data = datasets.make_fullranks_table()

    def test_make_fullranks_table(self):
        self.md5sum = joblib.hash(self.data)
        # print(self.md5sum)
        # assert self.md5sum == CHECKSUMS['fullranks']
        assert self.md5sum == '9af03c1ae884a9f663975bca2498e5cb'

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
        assert self.has_na == False, 'Fullranks table should not contain NA scores.' 

# END.
