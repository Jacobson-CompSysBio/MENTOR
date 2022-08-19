#!/usr/bin/env python
# coding: utf-8

'''
[TODO] Docstring
'''

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn import metrics
from scipy.spatial import distance
from scipy.cluster import hierarchy


def _root_mean_squared_error(y_true, y_pred=None, **kwargs):
    '''
    y_true, y_pred, *, sample_weight=None, multioutput='uniform_average', squared=True
    '''
    if y_pred is None:
        y_pred = np.linspace(y_true[0], y_true[-1], len(y_true))
    if len(y_true) != len(y_pred):
        raise ValueError('y_true and y_pred must be the same length.')
    mse = metrics.mean_squared_error(y_true, y_pred, **kwargs)
    rmse = np.sqrt(mse)
    return rmse


def _root_mean_squared_error_at_c(y_true, c, debug=False):
    b = len(y_true)
    l_start = 0
    l_stop = c+1
    r_start = c
    r_stop = b+1
    Lc = y_true[l_start:l_stop]
    Rc = y_true[r_start:r_stop]
    if debug:
        # print(l_start, l_stop, c, r_start, r_stop)
        # print(len(Lc), len(Rc))
        print(f'c={c}; Lc=[{l_start}:{l_stop}] ({len(Lc)}); Rc=[{r_start}:{r_stop}] ({len(Rc)})')
    rmse_c = ( ((c-1)/(b-1)) * _root_mean_squared_error(Lc) ) +              ( ((b-c)/(b-1)) * _root_mean_squared_error(Rc) )
    return rmse_c


def get_elbow(y_true, min_size=3, debug=False):
    '''
    RMSE_{c}={c-1\over b-1}\times RMSE(L_{c})+{b-c\over b-1}\times RMSE(R_{c}) \eqno{\hbox{[1]}}
    '''

    if isinstance(y_true, pd.DataFrame):
        raise ValueError('y_true must be a numpy array or pandas Series.')
    elif isinstance(y_true, pd.Series):
        y_true = y_true.values
    else:
        pass

    b = len(y_true)

    rmse_over_c = []

    for c in range(min_size, b-(min_size+1)):
        rmse_at_c = _root_mean_squared_error_at_c(y_true, c, debug=debug)
        rmse_over_c.append(rmse_at_c)
    # Adjust index by min_size.
    idx_of_elbow = int(np.argmin(rmse_over_c) + min_size)
    return idx_of_elbow


def cluster_hierarchical(path_or_dataframe):
    '''
    Parameters
    ----------
    path_or_dataframe : str, pd.DataFrame
        Path to 'fullranks' file from `RWR-CV --method=singletons` or
        pandas.DataFrame.

    Returns
    -------
    linkage_matrix
    '''
    if isinstance(path_or_dataframe, pd.DataFrame):
        fullranks = path_or_dataframe
    else:
        # Load the full ranks.
        fullranks = pd.read_table(path_or_dataframe)

    # Pivot full ranks -> ranks matrix.
    ranks = fullranks.pivot(index='seed', columns='NodeNames', values='rank')
    labels = ranks.index.to_list()

    # Find elbow and set max_rank.
    mean_scores = fullranks.groupby('rank')['Score'].mean()
    max_rank = get_elbow(mean_scores)

    # Filter the rank vectors.
    mask = (ranks <= max_rank).fillna(False)
    col_mask = mask.any()

    # Get pairwise distances.
    dmat = 1 - ranks.loc[:, col_mask].T.corr(method='spearman')
    dvec = distance.squareform(dmat)

    # Do the hierarchical clustering.
    Z = hierarchy.linkage(dvec, method='average')

    return Z, labels


def plot_dendrogram(Z, out_path=None, figsize='auto', draw_threshold=True, **kwargs):
    # Plot the dendrogram.
    if figsize == 'auto':
        width = 5
        height = np.shape(Z)[0] * 0.2
        if height < 10:
            height = 10
        figsize = (width, height)

    # if color_threshold == 'auto':
    #     kwargs['color_threshold'] = np.mean(Z[:,2])
    # elif color_threshold is None:
    #     # Default for hierarchy.dendrogram is to use `0.7*max(Z[:,2])`
    #     # when `color_threshold==None`. Set this explicitely to enable
    #     # `draw_threshold`, below.
    #     kwargs['color_threshold'] = 0.7*max(Z[:,2])
    # else:
    #     kwargs['color_threshold'] = float(color_threshold)

    plt.rc('figure', facecolor='white')
    plt.figure(figsize=figsize)

    if draw_threshold and kwargs.get('color_threshold') > 0:
        # You have to know the orientation: left/right > vline, top/bottom > hline.
        _orientation = kwargs.get('orientation', 'left')
        if _orientation == 'left' or _orientation == 'right':
            plot_line = plt.axvline
        elif _orientation == 'top' or _orientation == 'bottom':
            plot_line = plt.axhline
        else:
            raise ValueError(f'`orientation` must be one of ["top", "bottom", "left", "right"]: {_orientation}')

        plot_line(kwargs.get('color_threshold'), c='k', linewidth=1, linestyle='dotted')

    # One of the default colors for coloring the leaves is 'gray' (tab10 colors?).
    tree = hierarchy.dendrogram(
        Z,
        orientation=kwargs.pop('orientation', 'left'),
        labels=kwargs.pop('labels', None),
        leaf_font_size=kwargs.pop('leaf_font_size', 10),
        above_threshold_color=kwargs.pop('above_threshold_color', 'k'),
        count_sort=kwargs.pop('count_sort', True),
        **kwargs
    )

    if out_path is not None:
        plt.savefig(out_path, dpi=300, bbox_inches='tight')

    return tree


def get_clusters(Z, labels, threshold=0, match_to_leaves=None, out_path=None):
    '''
    Get clusters from Z at given threshold.

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
    clusters = hierarchy.fcluster(Z, t=threshold, criterion='distance')
    clusters = pd.DataFrame(zip(labels, clusters), columns=['label', 'cluster'])
    if match_to_leaves is not None:
        # The default dendrogram orientation is 'left'. Y-axis is 0 at the bottom, increasing.
        # So numeric tree['leaves'] start at 0 at the bottom of the figure.
        # To match the cluster labels to the dendrogram visually, reverse the
        # order of `tree['leaves']`.
        clusters = clusters.iloc[match_to_leaves[::-1]].reset_index(drop=True)

    if out_path is not None:
        clusters.to_csv(out_path, sep='\t', index=None)

    return clusters


def partition_fullranks(path_or_dataframe, out_dendrogram=None, out_clusters=None, threshold=0, **kwargs):
    Z, labels = cluster_hierarchical(path_or_dataframe)
    # if 'labels' not in kwargs:
    #     kwargs['labels'] = labels
    labels = kwargs.pop('labels', labels)

    if threshold == 'auto':
        threshold = np.mean(Z[:,2])
    elif threshold is None:
        # Default for hierarchy.dendrogram is to use `0.7*max(Z[:,2])` when
        # `color_threshold==None`. Set this explicitely to enable
        # `draw_threshold`, below.
        threshold = 0.7*max(Z[:,2])
    else:
        threshold = float(threshold)

    tree = plot_dendrogram(
        Z,
        labels=labels,
        color_threshold=threshold,
        out_path=out_dendrogram,
        **kwargs
    )
    clusters = get_clusters(
        Z,
        labels=labels,
        threshold=threshold,
        match_to_leaves=tree['leaves'],
        out_path=out_clusters
    )
    return tree, clusters


def parse_args(test=None):
    parser = argparse.ArgumentParser(
        description='Partition seeds from `RWR-CV --method=singletons ...` into clusters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # parser.add_argument(
    #     'positional',
    #     help='The positional argument.'
    # )
    parser.add_argument(
        '--rwr-fullranks', '-f',
        action='store',
        help='Path to "fullranks" file from `RWR-CV --method=singletons ...`'
    )
    parser.add_argument(
        '--partition', '-p',
        action='store_true',
        default=True,
        help='[PLACEHOLDER] Perform functional partitioning on "seed genes" from RWR fullranks file. This is the default.'
    )
    parser.add_argument(
        '--threshold', '-t',
        action='store',
        default=False,
        help='Perform functional partitioning on "seed genes" from RWR fullranks file.'
    )
    parser.add_argument(
        '--out-dendrogram', '-d',
        action='store',
        help='Save dendrogram to path.'
    )
    parser.add_argument(
        '--out-clusters', '-c',
        action='store',
        help='[TODO] NOT IMPLEMENTED. Save clusters to path.'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='count',
        default=0,
        help='Default: WARNING; once: INFO; twice: DEBUG'
    )

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if test is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(test)

    if args.verbose>0:
        # print(f"Arguments count: {len(args)}")
        for name, value in args.__dict__.items():
            print(f"{name}: {value}")

    return args


def main():
    args = parse_args()

    # Load the fullranks data once.
    fullranks = pd.read_table(args.rwr_fullranks)

    if args.partition:
        tree, clusters = partition_fullranks(
            path_or_dataframe=fullranks,
            threshold=args.threshold,
            out_dendrogram=args.out_dendrogram,
            out_clusters=args.out_clusters
        )

    return 0


if __name__ == '__main__':
    sys.exit(main())
