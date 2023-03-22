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

# Test linkage matrices. The linkage matrices returned from
# `scipy.hierarchy.linkage` and `HierarchicalClustering` should be the same.

class TestScipyClusteringRanks(unittest.TestCase):
    def setUp(self):
        _ranks = datasets.make_ranks_matrix()
        self.ranks = _ranks.fillna(0)
        self.dvec = distance.pdist(self.ranks, metric='euclidean')
        self.linkage_matrix = hierarchy.linkage(self.dvec, method='average')
        self.checksums = {
            'linkage_matrix': '4fd59486c7e1e78f2941c816b701ef61',
        }

    def hash_linkage_matrix(self):
        md5sum = joblib.hash(self.linkage_matrix)
        return md5sum

    def test_linkage_matrix(self):
        md5sum = self.hash_linkage_matrix()
        assert md5sum == self.checksums['linkage_matrix']


class TestScipyClusteringScores(unittest.TestCase):
    def setUp(self):
        _scores = datasets.make_scores_matrix()
        self.scores = _scores.fillna(_scores.max(axis=1))
        self.dvec = distance.pdist(self.scores, metric='euclidean')
        self.linkage_matrix = hierarchy.linkage(self.dvec, method='average')
        self.checksums = {
            'linkage_matrix': 'b1dd3e10cc1029f1a57202f4be3c9cfa',
        }

    def hash_linkage_matrix(self):
        md5sum = joblib.hash(self.linkage_matrix)
        return md5sum

    def test_linkage_matrix(self):
        md5sum = self.hash_linkage_matrix()
        assert md5sum == self.checksums['linkage_matrix']


class TestHierarchicalClusteringScoresEuclidean(unittest.TestCase):
    def setUp(self):
        # Load test data
        _scores = datasets.make_scores_matrix()
        self.scores = _scores.fillna(_scores.max(axis=1))
        self.mod = cluster.HierarchicalClustering(metric='euclidean')
        self.mod.fit(self.scores)
        self.checksums = {
            "clusters": "fb9a3f53b7eeb7c1f6e1db86108f8b5e",
            "linkage_matrix": "4817dbec688151d6a9956fb5ea936bd8"
        }

    def hash_linkage_matrix(self):
        md5sum = joblib.hash(self.mod.linkage_matrix)
        return md5sum

    def hash_clusters(self):
        md5sum = joblib.hash(self.mod.labels_)
        return md5sum

    def test_linkage_matrix(self):
        md5sum = self.hash_linkage_matrix()
        assert md5sum == self.checksums['linkage_matrix']

    def test_clusters(self):
        md5sum = self.hash_clusters()
        assert md5sum == self.checksums['clusters']


class TestHierarchicalClusteringScoresSpearman(unittest.TestCase):
    def setUp(self):
        # Load test data
        _scores = datasets.make_scores_matrix()
        self.scores = _scores.fillna(_scores.max(axis=1))
        self.mod = cluster.HierarchicalClustering(metric=metrics.spearman_d)
        self.mod.fit(self.scores)
        self.checksums = {
            "clusters": "fb9a3f53b7eeb7c1f6e1db86108f8b5e",
            "linkage_matrix": "bde826a41f0f9779eda84c340a68ee61"
        }

    def hash_linkage_matrix(self):
        md5sum = joblib.hash(self.mod.linkage_matrix)
        return md5sum

    def hash_clusters(self):
        md5sum = joblib.hash(self.mod.labels_)
        return md5sum

    def test_linkage_matrix(self):
        md5sum = self.hash_linkage_matrix()
        assert md5sum == self.checksums['linkage_matrix']

    def test_clusters(self):
        md5sum = self.hash_clusters()
        assert md5sum == self.checksums['clusters']


class TestHierarchicalClusteringRanksSpearman(unittest.TestCase):
    def setUp(self):
        # Load test data
        _ranks = datasets.make_ranks_matrix()
        self.ranks = _ranks.fillna(0)
        self.mod = cluster.HierarchicalClustering(metric=metrics.spearman_d)
        self.mod.fit(self.ranks)
        self.checksums = {
            "clusters": "fb9a3f53b7eeb7c1f6e1db86108f8b5e",
            "linkage_matrix": "f6d35fe51db2d93893c69d26f9535c3c"
        }

    def hash_linkage_matrix(self):
        md5sum = joblib.hash(self.mod.linkage_matrix)
        return md5sum

    def hash_clusters(self):
        md5sum = joblib.hash(self.mod.labels_)
        return md5sum

    def test_linkage_matrix(self):
        md5sum = self.hash_linkage_matrix()
        assert md5sum == self.checksums['linkage_matrix']

    def test_clusters(self):
        md5sum = self.hash_clusters()
        assert md5sum == self.checksums['clusters']

# END.
