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

# Test RWRtoolkit wrapper functions. Minimal testing here; just make sure the
# fullranks table is properly converted.

# [TODO] Add new test for elbow using scores-vs-ranks vector.

class TestRWRtoolkitFullranksTable(unittest.TestCase):
    def setUp(self):
        # Load test data
        self.fullranks = datasets.make_fullranks_table()
        self.scores_matrix = rwrtoolkit.fullranks_to_matrix(self.fullranks, to='scores')
        self.ranks_matrix = rwrtoolkit.fullranks_to_matrix(self.fullranks, to='ranks')
        self.checksums = {
            'fullranks_to_matrix_ranks': '88be7ef74f06a57127755af686ef312c',
            'fullranks_to_matrix_scores': '3c6d835ec34107f53bd13f9f1e1a16ce'
        }

    def hash_fullranks_to_matrix_scores(self):
        md5sum = joblib.hash(self.scores_matrix)
        return md5sum

    def hash_fullranks_to_matrix_ranks(self):
        md5sum = joblib.hash(self.ranks_matrix)
        return md5sum

    def test_fullranks_to_matrix_scores(self):
        md5sum = self.hash_fullranks_to_matrix_scores()
        assert md5sum == self.checksums['fullranks_to_matrix_scores']

    def test_fullranks_to_matrix_ranks(self):
        md5sum = self.hash_fullranks_to_matrix_ranks()
        assert md5sum == self.checksums['fullranks_to_matrix_ranks']

# END.
