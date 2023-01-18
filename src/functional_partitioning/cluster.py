'''
Cluster module

This file contains the HierarchicalClustering class, which is used to perform
hierarchical clustering on a set of data points.

The class is initialized with a matrix (samples in rows, features in columns)
and a distance function. The distance function should take two vectors as
input and return a scalar distance between them.
'''

import numpy as np
import pandas as pd
import logging

from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.utils.validation import check_memory
from sklearn.cluster import (
    AgglomerativeClustering,
    KMeans
)
from sklearn.cluster._agglomerative import _TREE_BUILDERS, _hc_cut
from sklearn.utils.validation import check_is_fitted
from scipy.cluster import hierarchy
from scipy.spatial import distance

LOGGER = logging.getLogger(__name__)

class HierarchicalClustering(AgglomerativeClustering, ClusterMixin, BaseEstimator):
    # Notes:
    # - `affinity` is passed to `pdist` as `metric`
    # - `affinity` can be a callable, but it takes a single matrix `X` as
    #   input, whereas `pdist` requires a callable takes two vectors `u` and
    #   `v` as arguments.

    def __init__(self, compute_linkage_matrix=True, compute_dendrogram=True, **kwargs):
        kwargs.setdefault('distance_threshold', 0)
        kwargs.setdefault('n_clusters', None)
        kwargs.setdefault('linkage', 'average')
        super().__init__(**kwargs)
        self.compute_linkage_matrix = compute_linkage_matrix
        self.compute_dendrogram = compute_dendrogram

    def _fit(self, X):
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
        # Copied from AgglomerativeClustering._fit >>>
        memory = check_memory(self.memory)

        if self.n_clusters is not None and self.n_clusters <= 0:
            raise ValueError(
                "n_clusters should be an integer greater than 0. %s was provided."
                % str(self.n_clusters)
            )

        if not ((self.n_clusters is None) ^ (self.distance_threshold is None)):
            raise ValueError(
                "Exactly one of n_clusters and "
                "distance_threshold has to be set, and the other "
                "needs to be None."
            )

        if self.distance_threshold is not None and not self.compute_full_tree:
            raise ValueError(
                "compute_full_tree must be True if distance_threshold is set."
            )

        if self.linkage == "ward" and self.affinity != "euclidean":
            raise ValueError(
                "%s was provided as affinity. Ward can only "
                "work with euclidean distances." % (self.affinity,)
            )

        if self.linkage not in _TREE_BUILDERS:
            raise ValueError(
                "Unknown linkage type %s. Valid options are %s"
                % (self.linkage, _TREE_BUILDERS.keys())
            )
        tree_builder = _TREE_BUILDERS[self.linkage]

        connectivity = self.connectivity
        if self.connectivity is not None:
            if callable(self.connectivity):
                connectivity = self.connectivity(X)
            connectivity = check_array(
                connectivity, accept_sparse=["csr", "coo", "lil"]
            )

        n_samples = len(X)
        compute_full_tree = self.compute_full_tree
        if self.connectivity is None:
            compute_full_tree = True
        if compute_full_tree == "auto":
            if self.distance_threshold is not None:
                compute_full_tree = True
            else:
                # Early stopping is likely to give a speed up only for
                # a large number of clusters. The actual threshold
                # implemented here is heuristic
                compute_full_tree = self.n_clusters < max(100, 0.02 * n_samples)
        n_clusters = self.n_clusters
        if compute_full_tree:
            n_clusters = None

        # Construct the tree
        kwargs = {}
        if self.linkage != "ward":
            kwargs["linkage"] = self.linkage
            kwargs["affinity"] = self.affinity

        distance_threshold = self.distance_threshold

        return_distance = (distance_threshold is not None) or self.compute_distances

        out = memory.cache(tree_builder)(
            X,
            connectivity=connectivity,
            n_clusters=n_clusters,
            return_distance=return_distance,
            **kwargs,
        )
        (self.children_, self.n_connected_components_, self.n_leaves_, parents) = out[
            :4
        ]

        if return_distance:
            self.distances_ = out[-1]

        if self.distance_threshold is not None:  # distance_threshold is used
            self.n_clusters_ = (
                np.count_nonzero(self.distances_ >= distance_threshold) + 1
            )
        else:  # n_clusters is used
            self.n_clusters_ = self.n_clusters

        # Cut the tree
        if compute_full_tree:
            self.labels_ = _hc_cut(self.n_clusters_, self.children_, self.n_leaves_)
        else:
            labels = _hierarchical.hc_get_heads(parents, copy=False)
            # copy to avoid holding a reference on the original array
            labels = np.copy(labels[:n_samples])
            # Reassign cluster numbers
            self.labels_ = np.searchsorted(np.unique(labels), labels)
        # End of AgglomerativeClustering._fit <<<
            
        # Compute the linkage matrix (**NEW**).
        if self.compute_linkage_matrix or self.compute_dendrogram:
            # Create the counts of samples under each node.
            # https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html#sphx-glr-auto-examples-cluster-plot-agglomerative-dendrogram-py
            counts = np.zeros(self.children_.shape[0])
            n_samples = len(self.labels_)
            for i, merge in enumerate(self.children_):
                current_count = 0
                for child_idx in merge:
                    if child_idx < n_samples:
                        current_count += 1  # leaf node
                    else:
                        current_count += counts[child_idx - n_samples]
                counts[i] = current_count

            self.linkage_matrix = np.column_stack(
                [self.children_, self.distances_, counts]
            ).astype(float)
        else:
            self.linkage_matrx = None
            
        # Compute the dendrogram (**NEW**).
        if self.compute_dendrogram:
            self.dendrogram = hierarchy.dendrogram(self.linkage_matrix, no_plot=True)
        else:
            self.dendrogram = None
    

        return self


########################################################################
# DEPRECATED: The functions below will be removed in a future version.
########################################################################


def cluster_hierarchical(X, corr_method='spearman', linkage_method='average'):
    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X)

    # Get pairwise distances.
    dmat = 1 - X.T.corr(method=corr_method)
    dvec = distance.squareform(dmat)

    # Do the hierarchical clustering.
    Z = hierarchy.linkage(dvec, method=linkage_method)

    # return Z, labels
    return Z


def get_clusters(Z, labels=None, threshold=None, n_clusters=None, match_to_leaves=None, out_path=None):
    '''
    Get clusters from linkage matrix `Z` at given threshold.

    Parameters
    ----------
    Z : linkage matrix
    labels : list
    threshold : int, float
    match_to_leave : None, list
        List of ints to order the rows. Use `tree['leaves']` to visually match
        the list of cluster labels to the leaves of the dendrogram.  See
        example below.

    Examples
    --------
    Z, labels = cluster_hierarchical(fullranks)
    tree = plot_dendrogram(Z, labels)

    # Genes are ordered in same order as `labels`:
    clusters = get_clusters(Z, labels, threshold=t)

    # Genes are ordered the same as `tree['leaves']`:
    clusters = get_clusters(Z, labels, threshold=t, match_to_leaves=tree['leaves'])
    '''
    # clusters = hierarchy.fcluster(Z, t=threshold, criterion='distance')
    # clusters = pd.DataFrame(zip(labels, clusters), columns=['label', 'cluster'])
    
    # There is a discrepancy bn `fcluster` and `cut_tree`, when the cut is on the branch,
    # one returns the clustering *above* the branch but the other returns the clustering *below*
    # the branch (don't recall off hand which is which).
    # The two methods (`fcluster` vs `cut_tree`) differ in implementation:
    # >>> threshold = 0.3
    # >>> hierarchy.fcluster(Z, t=threshold, criterion='distance')
    # array([2, 2, 1, 2, 1, 1, 1, 3, 1], dtype=int32)
    # >>> hierarchy.cut_tree(Z, height=threshold, n_clusters=None)
    # array([[0],
    #        [0],
    #        [1],
    #        [0],
    #        [1],
    #        [1],
    #        [1],
    #        [2],
    #        [1]])

    if not threshold and not n_clusters:
        # Get a matrix of all clusterings, one for each agglomeration step.
        clusters = hierarchy.cut_tree(Z, n_clusters=None, height=None)
    else:
        # Get a vector of a single clustering.
        clusters = hierarchy.cut_tree(Z, n_clusters=n_clusters, height=threshold).flatten()
    clusters = pd.DataFrame(clusters, index=labels)
    clusters.index.name = 'node_label'

    if match_to_leaves is not None:
        # The default dendrogram orientation is 'left'. Y-axis is 0 at the bottom, increasing.
        # So numeric tree['leaves'] start at 0 at the bottom of the figure.
        # To match the cluster labels to the dendrogram visually, reverse the
        # order of `tree['leaves']`.
        # clusters = clusters.iloc[match_to_leaves[::-1]].reset_index(drop=True)
        clusters = clusters.iloc[match_to_leaves[::-1]]

    if out_path is not None:
        clusters.to_csv(out_path, sep='\t')
        LOGGER.info(f'Saved clusters: {out_path}')

    return clusters

# END
