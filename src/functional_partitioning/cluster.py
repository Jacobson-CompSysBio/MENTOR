'''
Cluster module

This file contains the HierarchicalClustering class, which is used to perform
hierarchical clustering on a set of data points.

The class is initialized with a matrix (samples in rows, features in columns)
and a distance function. The distance function should take two vectors as
input and return a scalar distance between them.
'''

import dynamicTreeCut
import logging
import numpy as np
import pandas as pd
import warnings

from dynamicTreeCut import R_func as dtc_utils
from scipy.cluster import hierarchy
from scipy.spatial import distance
from sklearn import metrics
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.cluster import ( AgglomerativeClustering, KMeans )
from sklearn.cluster._agglomerative import _TREE_BUILDERS, _hc_cut
from sklearn.utils.validation import check_is_fitted, check_memory


LOGGER = logging.getLogger(__name__)


def check_symmetry(dmat, atol=1e-6):
    # Manually check the distance matrix for symmetry with absolute value tolerance.
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
        n_clusters=None, # Changed default.
        *,
        metric='euclidean',
        memory=None,
        connectivity=None,
        compute_full_tree="auto",
        linkage="average",
        distance_threshold=None,
        compute_distances=False,
        compute_linkage_matrix=True,
        compute_dendrogram=True,
        cut_threshold=None,
        cut_method='cutreeHybrid',
    ):
        # # Handle 'affinity', which will be removed in v1.4. Sklearn passes
        # # 'affinity' to 'pdist' as 'metric' (eg, in the `_average` method).
        # # After v1.4, I assume this will change so that 'metric' will passed
        # # instead of 'affinity'.
        # if metric is None:
        #     if affinity is None or affinity == 'deprecated':
        #         metric = 'euclidean'
        #     else:
        #         metric = affinity
        super_kwargs = dict(
            affinity = metric,
            metric = metric,
            n_clusters = n_clusters,
            distance_threshold = distance_threshold,
            memory = memory,
            connectivity = connectivity,
            compute_full_tree = compute_full_tree,
            linkage = linkage,
            compute_distances = compute_distances,
        )
        super().__init__(**super_kwargs)
        # [IN_PROGRESS] Init all the attributes ... 
        self._cutreeHybrid = None
        self._n_samples = None
        self.compute_dendrogram = compute_dendrogram
        self.compute_linkage_matrix = compute_linkage_matrix
        self.cut_method = cut_method
        self.cut_threshold = cut_threshold
        self.dendrogram = None
        self.distances_ = 'deprecated'
        self.labels_ = None  # Not sure if this is needed.
        self.method = linkage # [TODO] Use method instead of linkage.
        self.min_cluster_size_ = None
        # [TODO] Change 'linkage' to 'linkage_method'
        self.linkage_method = None
        # [TODO] Change 'metric' to 'linkage_metric'
        self.linkage_metric = None

    # [TODO] Add fit method
    # [TODO] Add fit_predict method
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
        # [TODO] What does `check_memory` do?
        memory = check_memory(self.memory)

        if self.cut_method != 'cutreeHybrid':
            # Only one of n_clusters or distance_threshold can be passed to hierarchy.cut_tree.
            if self.n_clusters is not None and self.distance_threshold is not None:
                raise ValueError("Provide n_clusters OR distance_threshold, not both.")
            if self.n_clusters is not None and self.n_clusters <= 0:
                raise ValueError(
                    "n_clusters should be an integer greater than 0. %s was provided."
                    % str(self.n_clusters)
                )

        # Check for consistency between shape of X and the distance metric.
        if X.ndim == 1 and self.metric != "precomputed":
            raise ValueError("X should be a 2D array if metric is \"%s\"." % self.metric)
        # # Removed: I think scikit learn only wants a 2D array.
        # elif X.ndim != 1 and self.metric == "precomputed":
        #     raise ValueError("X should be a 1D condensed distance matrix if metric is \"precomputed\"." % self.metric)

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

        # Compute pairwise distances between samples.
        if X.ndim == 1 and self.metric == 'precomputed':
            # It's a condensed distance matrix (ie, vector), but scikit-learn
            # pairwise_distances function wants it as a square, uncondensed
            # distance matrix.
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

        # # Manually check the distance matrix for symmetry with absolute value tolerance.
        # atol = 1e-6
        # assert dmat.shape[0] == dmat.shape[1], f'The distance matrix is not square: {dmat.shape}'
        # np.testing.assert_allclose(
        #     np.diag(dmat),
        #     0,
        #     atol=atol,
        #     err_msg=f'The diagonal values are not within tolerance to 0 (`diag - 0 > {atol}`)'
        # )
        # np.testing.assert_allclose(
        #     dmat - dmat.T,
        #     0,
        #     atol=atol,
        #     err_msg=f'The distance matrix is not symmetric (`d - d.T >{atol}`)'
        # )
        is_symmetric = check_symmetry(dmat)
        self.pairwise_distances = distance.squareform(dmat, checks=False)

        # Compute the linkage matrix.
        self.linkage_matrix = hierarchy.linkage(
            self.pairwise_distances,
            metric='precomputed',
            method=self.linkage,
        )
        
        # Get the clusters.
        if self.cut_method == 'cutreeHybrid':
            # self.min_cluster_size_ = int(np.sqrt(self._n_samples))
            self.min_cluster_size_ = 3
            # Set the cutHeight parameter using default method (99% of dendrogram height).
            # Passing this parameter is necessary to avoid a printed message from dynamicTreeCut.
            # https://github.com/kylessmith/dynamicTreeCut/blob/3734243ee547bb9c220e5aef046587ca1694c7a7/dynamicTreeCut/dynamicTreeCut.py#L176
            dendro_height = dtc_utils.get_heights(self.linkage_matrix)
            dendro_merge = dtc_utils.get_merges(self.linkage_matrix)
            nMerge = len(dendro_height)
            refQuantile = 0.05
            refMerge = np.round(nMerge * refQuantile)
            if refMerge < 1:
                refMerge = 1
            refHeight = dendro_height[int(refMerge) - 1]
            self.cutHeight_ = 0.99 * (np.max(dendro_height) - refHeight) + refHeight

            self._cutreeHybrid = dynamicTreeCut.cutreeHybrid(
                self.linkage_matrix,
                self.pairwise_distances,
                minClusterSize=self.min_cluster_size_,
                cutHeight=self.cutHeight_,
                verbose=0,
                # Defaults:
                deepSplit=1,
                maxCoreScatter=None,
                minGap=None,
                maxAbsCoreScatter=None,
                minAbsGap=None,
                minSplitHeight=None,
                minAbsSplitHeight=None,
                externalBranchSplitFnc=None,
                minExternalSplit=None,
                externalSplitOptions=[],
                externalSplitFncNeedsDistance=None,
                assumeSimpleExternalSpecification=True,
                pamStage=True,
                pamRespectsDendro=True,
                useMedoids=False,
                maxPamDist=None,
                respectSmallClusters=True,
                indent=0,
            )
            self.labels_ = self._cutreeHybrid['labels']
        else:
            self.labels_ = hierarchy.cut_tree(
                self.linkage_matrix,
                n_clusters=self.n_clusters,
                height=self.distance_threshold
            )

        if self.labels_.ndim == 1:
            # This really only applies to cutreeHybrid.
            self.labels_ = self.labels_.reshape(-1, 1)

        # Compute the dendrogram.
        if self.compute_dendrogram:
            self.dendrogram = hierarchy.dendrogram(self.linkage_matrix, no_plot=True)
    
        return self


########################################################################
# DEPRECATED: The functions below will be removed in a future version. {{{
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
# }}}
########################################################################
# END
