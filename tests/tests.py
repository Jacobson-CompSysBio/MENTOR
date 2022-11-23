#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import string
import joblib

# from sklearn.datasets import make_blobs
from sklearn import datasets, preprocessing
from matplotlib import pyplot


from functional_partitioning import functional_partitioning as fp


CHECKSUMS = {
    'features': '05f5c7599ac919bddb182aec8d2aaf16',
    'scores_matrix': '0bc8e7374e6aa2b56ca100eabad48e4b',
    'ranks_matrix': '946e98267a27f4d837bdf0cbb4d1a529',
    'fullranks': 'e8646ece5873b493e38a71f877372e66'
}


RANDOM_STATE = 42


def make_features(random_state=42):
    scaler = preprocessing.MinMaxScaler(feature_range=(0.01, 0.99))
    X, y = datasets.make_blobs(n_samples=9, n_features=26, centers=2, random_state=random_state)
    # print(X.min(), X.max())
    X = scaler.fit_transform(X)
    # print(X.min(), X.max())
    for i in range(X.shape[0]):
        X[i, i] = np.nan
    # print(X.min(), X.max())
    X.shape
    return X, y

def test_make_features():
    X, y = make_features()
    assert joblib.hash(X) == CHECKSUMS['features']


def make_scores_matrix(random_state=42):
    X, y = make_features()
    scores = pd.DataFrame(X, index=list(string.ascii_uppercase[:X.shape[0]]), columns=list(string.ascii_uppercase[:X.shape[1]]))
    return scores

def test_make_scores_matrix():
    scores_matrix = make_scores_matrix()
    assert joblib.hash(scores_matrix) == CHECKSUMS['scores_matrix']


def make_ranks_matrix(random_state=42):
    scores = make_scores_matrix()
    ranks = scores.rank(axis=1, ascending=False)
    return ranks

def test_make_ranks_matrix():
    ranks_matrix = make_ranks_matrix()
    assert joblib.hash(ranks_matrix) == CHECKSUMS['ranks_matrix']


def make_fullranks_table(random_state=42):
    scores = make_scores_matrix(random_state=random_state)
    ranks = make_ranks_matrix(random_state=random_state)
    fullranks = pd.merge(
        left=scores.stack().to_frame(name='Score'),
        right=ranks.stack().to_frame(name='rank'),
        left_index=True,
        right_index=True
    ).reset_index().rename(columns={'level_0': 'seed', 'level_1': 'NodeNames'}).dropna(how='any')
    return fullranks

def test_make_fullranks_table():
    fullranks = make_fullranks_table()
    assert joblib.hash(fullranks) == CHECKSUMS['fullranks']

# END.
