'''
Cluster module

This file contains the HierarchicalClustering class, which is used to perform
hierarchical clustering on a set of data points.

The class is initialized with a matrix (samples in rows, features in columns)
and a distance function. The distance function should take two vectors as
input and return a scalar distance between them.
'''

import logging
import numpy as np
import pandas as pd
import warnings
from dynamicTreeCut import dynamicTreeCut
from dynamicTreeCut import R_func as dtc_utils
from scipy.cluster import hierarchy
from scipy.spatial import distance
from mentor import _metrics as metrics
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.cluster import ( AgglomerativeClustering, KMeans )
from sklearn.cluster._agglomerative import _TREE_BUILDERS, _hc_cut
from sklearn.utils.validation import check_is_fitted, check_memory

LOGGER = logging.getLogger(__name__)

def check_symmetry(dmat, atol=1e-6):
    assert dmat.shape[0] == dmat.shape[1], f'The distance matrix is not square: {dmat.shape}'
    np.testing.assert_allclose(
        np.diag(dmat),
        0,
        atol=atol,
        err_msg=f'The diagonal values are not within tolerance to 0 (`diag - 0 > {atol}`)'
    )
    np.testing.assert_allclose(
        dmat - dmat.T,
        0,
        atol=atol,
        err_msg=f'The distance matrix is not symmetric (`d - d.T >{atol}`)'
    )
    return True

class HierarchicalClustering(AgglomerativeClustering):
    # Notes:
    # - `affinity` is deprecated in sklearn v1.2 and will be removed in v1.4.
    # - `affinity` is passed to `pdist` as `metric`
    # - `affinity` can be a callable, but it takes a single matrix `X` as
    #   input, whereas `pdist` requires a callable takes two vectors `u` and
    #   `v` as arguments.
    def __init__(
        self,
        *,
        metric='euclidean',
        memory=None,
        connectivity=None,
        linkage="average",
        compute_distances=False,
        compute_linkage_matrix=True
    ):
        super_kwargs = dict(
            affinity = metric,
            metric = metric,
            memory = memory,
            connectivity = connectivity,
            linkage = linkage,
            compute_distances = compute_distances,
        )
        super().__init__(**super_kwargs)
        self._n_samples = None
        self.compute_linkage_matrix = compute_linkage_matrix
        self.distances_ = 'deprecated'
        self.method = linkage # use <method> instead of linkage
        self.linkage_method = None
        self.linkage_metric = None # use <linkage_metric> instead

    def _fit(self, X, **kwargs):
        """
        Fit without validation
        
        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features) or (n_samples, n_samples)
            Training instances to cluster, or distances between instances if
            ``affinity='precomputed'``.
        Returns
        -------
        self : object
            Returns the fitted instance.
            
        https://github.com/scikit-learn/scikit-learn/blob/f3f51f9b6/sklearn/cluster/_agglomerative.py#L917
        """
        memory = check_memory(self.memory)
        if X.ndim == 1 and self.metric != "precomputed":
            raise ValueError("X should be a 2D array if metric is \"%s\"." % self.metric)

        if self.linkage == 'ward' and self.metric == 'precomputed':
            warnings.warn(
                    'The Ward linkage algorithm is only valid for Euclidean distances; you gave `precomputed`'
                    ' Be sure that your distance metric is Euclidean.'
            )
        elif self.linkage == "ward" and self.metric != "euclidean":
            raise ValueError(
                'The Ward linkage algorithm is only valid for Euclidean distances; you gave `%s`.'
                % self.metric
            )
        if X.ndim == 1 and self.metric == 'precomputed':
            features = distance.squareform(X, checks=False)
        else:
            features = X
        self._n_samples = features.shape[0]
        pairwise_distances_kwargs = dict(
            Y=kwargs.get('Y', None),
            metric=self.metric,
            n_jobs=kwargs.get('n_jobs', None),
            force_all_finite=kwargs.get('force_all_finite', True)
        )
        dmat = metrics.pairwise_distances(features, **pairwise_distances_kwargs)
        is_symmetric = check_symmetry(dmat)
        self.pairwise_distances = distance.squareform(dmat, checks=False)
        #self.linkage_matrix = hierarchy.linkage(
        #    self.pairwise_distances,
        #    metric='precomputed',
        #    method=self.linkage,
        #    optimal_ordering=self.optimal_ordering
        #)
        return self
