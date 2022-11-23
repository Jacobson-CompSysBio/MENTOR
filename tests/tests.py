#!/usr/bin/env python
# coding: utf-8

import joblib

from sklearn import datasets, preprocessing
from matplotlib import pyplot

from functional_partitioning import functional_partitioning as fp
from functional_partitioning import datasets


CHECKSUMS = {
    'features': '05f5c7599ac919bddb182aec8d2aaf16',
    'scores_matrix': '0bc8e7374e6aa2b56ca100eabad48e4b',
    'ranks_matrix': '946e98267a27f4d837bdf0cbb4d1a529',
    'fullranks': 'e8646ece5873b493e38a71f877372e66'
}


RANDOM_STATE = 42


def test_make_features():
    X, y = datasets.make_features()
    assert joblib.hash(X) == CHECKSUMS['features']


def test_make_scores_matrix():
    scores_matrix = datasets.make_scores_matrix()
    assert joblib.hash(scores_matrix) == CHECKSUMS['scores_matrix']


def test_make_ranks_matrix():
    ranks_matrix = datasets.make_ranks_matrix()
    assert joblib.hash(ranks_matrix) == CHECKSUMS['ranks_matrix']


def test_make_fullranks_table():
    fullranks = datasets.make_fullranks_table()
    assert joblib.hash(fullranks) == CHECKSUMS['fullranks']


# END.
