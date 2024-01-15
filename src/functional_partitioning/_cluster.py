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
from functional_partitioning import _metrics as metrics
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


#def calc_cut_threshold(linkage_matrix, threshold, features=None):
#    '''
#    Calculate a threshold for `cut_tree` from the linkage matrix `linkage_matrix`.
#    '''
#    if threshold == 'mean':
#        threshold = np.mean(linkage_matrix[:,2])
#    elif threshold == 'best_chi':
#        if features is None:
#            raise ValueError('`features` must be provided if `threshold` is "best_chi"')
#        clusterings = hierarchy.cut_tree(linkage_matrix, n_clusters=None, height=None)
#        chi_features = metrics.calc_chi(features, clusterings)
#        best_at = np.nan_to_num(chi_features).argmax()
#        h1 = linkage_matrix[best_at, 2]
#        h0 = linkage_matrix[best_at-1, 2]
#        threshold = np.mean((h0, h1))
#    else:
#        pass

#    return threshold


class HierarchicalClustering(AgglomerativeClustering):
    # Notes:
    # - `affinity` is deprecated in sklearn v1.2 and will be removed in v1.4.
    # - `affinity` is passed to `pdist` as `metric`
    # - `affinity` can be a callable, but it takes a single matrix `X` as
    #   input, whereas `pdist` requires a callable takes two vectors `u` and
    #   `v` as arguments.

    def __init__(
        self,
        #n_clusters=None, # Changed default.
        *,
        metric='euclidean',
        memory=None,
        connectivity=None,
        #compute_full_tree="auto",
        linkage="average",
        #distance_threshold=None,
        compute_distances=False,
        compute_linkage_matrix=True
        #compute_dendrogram=True,
        #cut_threshold=None,
        #cut_method='cutreeHybrid',
    ):

        super_kwargs = dict(
            affinity = metric,
            metric = metric,
            #n_clusters = n_clusters,
            #distance_threshold = distance_threshold,
            memory = memory,
            connectivity = connectivity,
            #compute_full_tree = compute_full_tree,
            linkage = linkage,
            compute_distances = compute_distances,
        )
        super().__init__(**super_kwargs)
        #self._cutreeHybrid = None
        self._n_samples = None
        #self.compute_dendrogram = compute_dendrogram
        self.compute_linkage_matrix = compute_linkage_matrix
        #self.cut_method = cut_method
        #self.cut_threshold = cut_threshold
        #self.dendrogram = None
        self.distances_ = 'deprecated'
        #self.labels_ = None
        self.method = linkage # use <method> instead of linkage
        #self.min_cluster_size_ = None
        self.linkage_method = None
        self.linkage_metric = None # use <linkage_metric> instead

    # add <fit> and <fit_predict> method
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
        
        #if self.distance_threshold is not None:
        #    warnings.warn(
        #        "The parameter 'distance_threshold' is deprecated in 0.5.1 "
        #        "and will be removed in a future version. Use `cut_threshold`.",
        #        DeprecationWarning,
        #    )
        #    if self.cut_threshold is None:
        #        self.cut_threshold = self.distance_threshold

        #if self.cut_method != 'cutreeHybrid':
        #    if self.n_clusters is not None and self.cut_threshold is not None:
        #        raise ValueError("Provide n_clusters OR cut_threshold, not both.")
        #    if self.n_clusters is not None and self.n_clusters <= 0:
        #        raise ValueError(
        #            "n_clusters should be an integer greater than 0. %s was provided."
        #            % str(self.n_clusters)
        #        )

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
        
        #if self.cut_method == 'cutreeHybrid':
        #    self.min_cluster_size_ = 1
        #    dendro_height = dtc_utils.get_heights(self.linkage_matrix)
        #    dendro_merge = dtc_utils.get_merges(self.linkage_matrix)
        #    nMerge = len(dendro_height)
        #    refQuantile = 0.05
        #    refMerge = np.round(nMerge * refQuantile)
        #    if refMerge < 1:
        #        refMerge = 1
        #    refHeight = dendro_height[int(refMerge) - 1]
        #    self.cutHeight_ = 0.99 * (np.max(dendro_height) - refHeight) + refHeight

            #self._cutreeHybrid = dynamicTreeCut.cutreeHybrid(
                #self.linkage_matrix,
                #self.pairwise_distances,
                #minClusterSize=self.min_cluster_size_,
                #cutHeight=self.cutHeight_,
                #verbose=0,
                #deepSplit=1,
                #maxCoreScatter=None,
                #minGap=None,
                #maxAbsCoreScatter=None,
                #minAbsGap=None,
                #minSplitHeight=None,
                #minAbsSplitHeight=None,
                #externalBranchSplitFnc=None,
                #minExternalSplit=None,
                #externalSplitOptions=[],
                #externalSplitFncNeedsDistance=None,
                #assumeSimpleExternalSpecification=True,
                #pamStage=True,
                #pamRespectsDendro=True,
                #useMedoids=False,
                #maxPamDist=None,
                #respectSmallClusters=True,
                #indent=0,
            #)
            #self.labels_ = self._cutreeHybrid['labels']
            #self.cut_threshold_ = None
        #else:
            #self.cut_threshold_ = calc_cut_threshold(
                #self.linkage_matrix,
                #self.cut_threshold,
                #features=features
            #)
            #self.labels_ = hierarchy.cut_tree(
                #self.linkage_matrix,
                #n_clusters=self.n_clusters,
                #height=self.cut_threshold_
            #)

        #if self.labels_.ndim == 1:
            #self.labels_ = self.labels_.reshape(-1, 1)

        #if self.compute_dendrogram:
            #self.dendrogram = hierarchy.dendrogram(self.linkage_matrix, no_plot=True)
    
        return self


########################################################################
# DEPRECATED: The functions below will be removed in a future version
########################################################################


#def cluster_hierarchical(X, corr_method='spearman', linkage_method='average'):
#    if not isinstance(X, pd.DataFrame):
#        X = pd.DataFrame(X)
#    dmat = 1 - X.T.corr(method=corr_method)
#    dvec = distance.squareform(dmat)
#    Z = hierarchy.linkage(dvec, method=linkage_method)
#    return Z


#def get_clusters(Z, labels=None, threshold=None, n_clusters=None, match_to_leaves=None, out_path=None):
#    '''
#    Get clusters from linkage matrix `Z` at given threshold.
#
#    Parameters
#    ----------
#    Z : linkage matrix
#    labels : list
#    threshold : int, float
#    match_to_leave : None, list
#        List of ints to order the rows. Use `tree['leaves']` to visually match
#        the list of cluster labels to the leaves of the dendrogram.  See
#        example below.
#
#    Examples
#    --------
#    Z, labels = cluster_hierarchical(fullranks)
#    tree = plot_dendrogram(Z, labels)
#
#    # Genes are ordered in same order as `labels`:
#    clusters = get_clusters(Z, labels, threshold=t)
#
#    # Genes are ordered the same as `tree['leaves']`:
#    clusters = get_clusters(Z, labels, threshold=t, match_to_leaves=tree['leaves'])
#    '''
#    if not threshold and not n_clusters:
#        clusters = hierarchy.cut_tree(Z, n_clusters=None, height=None)
#    else:
#        clusters = hierarchy.cut_tree(Z, n_clusters=n_clusters, height=threshold).flatten()
#    clusters = pd.DataFrame(clusters, index=labels)
#    clusters.index.name = 'node_label'
#
#    if match_to_leaves is not None:
#        clusters = clusters.iloc[match_to_leaves[::-1]]
#
#    if out_path is not None:
#        clusters.to_csv(out_path, sep='\t')
#        LOGGER.info(f'Saved clusters: {out_path}')
#
#    return clusters

########################################################################
