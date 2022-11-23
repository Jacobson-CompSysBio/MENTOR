import pandas as pd
import numpy as np
import string
import joblib

# from sklearn.datasets import make_blobs
from sklearn import datasets, preprocessing
from matplotlib import pyplot

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

def make_scores_matrix(random_state=42):
    X, y = make_features()
    scores = pd.DataFrame(X, index=list(string.ascii_uppercase[:X.shape[0]]), columns=list(string.ascii_uppercase[:X.shape[1]]))
    return scores

def make_ranks_matrix(random_state=42):
    scores = make_scores_matrix()
    ranks = scores.rank(axis=1, ascending=False)
    return ranks

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

