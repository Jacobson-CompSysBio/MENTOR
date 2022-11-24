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
import logging

from sklearn import metrics
from scipy.spatial import distance
from scipy.cluster import hierarchy

# from .version import __version__
# Error when run as a script:
# > ImportError: attempted relative import with no known parent package

DPI = 300
LOGGER = logging.getLogger(__name__)

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
    LOGGER.debug(f'c={c}; Lc=[{l_start}:{l_stop}] ({len(Lc)}); Rc=[{r_start}:{r_stop}] ({len(Rc)})')
    rmse_c = ( ((c-1)/(b-1)) * _root_mean_squared_error(Lc) ) +              ( ((b-c)/(b-1)) * _root_mean_squared_error(Rc) )
    return rmse_c


def get_elbow(y_true, min_size=3, debug=False):
    '''
    RMSE_{c}={c-1\over b-1}\times RMSE(L_{c})+{b-c\over b-1}\times RMSE(R_{c}) \eqno{\hbox{[1]}}

    Parameters
    ----------
    y_true : array-like
        The true values.
    min_size : int
        Minimum size of the left and right clusters.
    debug : bool
        Print debug information.

    Returns
    -------
    c : int
        The index of the elbow.

    Examples
    --------
    >>> y_true = scores_matrix.values
    >>> c = get_elbow(y_true)
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


def cluster_hierarchical(path_or_dataframe, max_rank='elbow', drop_missing=True):
    '''
    Parameters
    ----------
    path_or_dataframe : str, pd.DataFrame
        Path to 'fullranks' file from `RWR-CV --method=singletons` or
        pandas.DataFrame.
    max_rank : int, str
        Maximum rank to use for clustering. If 'elbow', use elbow method to
        determine max_rank.
    drop_missing : bool
        Drop genes that are labeled "missing" in the fullranks file.

    Returns
    -------
    linkage_matrix, labels
    '''
    if isinstance(path_or_dataframe, pd.DataFrame):
        fullranks = path_or_dataframe
    else:
        # Load the full ranks.
        fullranks = pd.read_table(path_or_dataframe)

    if drop_missing:
        # Drop rows with seeds that were not found in the network.
        fullranks = fullranks.query('seed!="missing"')

    # Pivot full ranks -> ranks matrix.
    ranks = fullranks.pivot(index='seed', columns='NodeNames', values='rank')
    labels = ranks.index.to_list()

    if max_rank == 'elbow':
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

    # return Z, labels
    return dict(linkage=Z, labels=labels, max_rank=max_rank)


def savefig(out_path=None, **kwargs):
    kwargs.setdefault('dpi', DPI)
    kwargs.setdefault('bbox_inches', 'tight')
    if out_path is not None:
        try:
            plt.savefig(out_path, **kwargs)
            LOGGER.info('Saved figure: %s', out_path)
        except Exception as e:
            LOGGER.error('Failed to save figure: %s', str(e))


def plot_dendrogram(
    Z,
    out_path=None,
    figsize='auto',
    draw_threshold=True,
    **kwargs
):
    '''
    labels : None
        Dummy (for consistency with `plot_dendrogram_polar`).
    out_path : str
    figsize : tuple, str
        Default is 'auto'.
    draw_threshold : bool
    kwargs
        Passed to `hierarchy.dendrogram`
    '''
    if kwargs.get('no_plot'):
        pass
    elif plt.get_fignums() and kwargs.get('ax') is None:
        # A figure exists; use it.
        # > To test whether there is currently a figure on the pyplot figure
        # > stack, check whether `~.pyplot.get_fignums()` is empty.
        # > ~ help(plt.gcf)
        kwargs['ax'] = plt.gca()
    elif kwargs.get('ax') is None:
        if figsize == 'auto':
            width = 5
            height = np.shape(Z)[0] * 0.2
            if height < 10:
                height = 10
            figsize = (width, height)
        # Initialize the figure.
        plt.rc('figure', facecolor='white')
        plt.figure(figsize=figsize)

    if kwargs.get('no_plot'):
        pass
    elif draw_threshold and kwargs.get('color_threshold', 0) > 0:
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
        p=kwargs.get('p', 30),
        truncate_mode=kwargs.get('truncate_mode', None),
        color_threshold=kwargs.get('color_threshold', None),
        get_leaves=kwargs.get('get_leaves', True),
        orientation=kwargs.get('orientation', 'left'),
        labels=kwargs.get('labels', None),
        count_sort=kwargs.get('count_sort', True),
        distance_sort=kwargs.get('distance_sort', False),
        show_leaf_counts=kwargs.get('show_leaf_counts', True),
        no_plot=kwargs.get('no_plot', False),
        no_labels=kwargs.get('no_labels', False),
        leaf_font_size=kwargs.get('leaf_font_size', 10),
        leaf_rotation=kwargs.get('leaf_rotation', None),
        leaf_label_func=kwargs.get('leaf_label_func', None),
        show_contracted=kwargs.get('show_contracted', False),
        link_color_func=kwargs.get('link_color_func', None),
        ax=kwargs.get('ax', None),
        above_threshold_color=kwargs.get('above_threshold_color', 'k')
    )

    if kwargs.get('no_plot'):
        pass
    else:
        savefig(out_path=out_path)

    return tree


def plot_dendrogram_polar(
    Z,
    labels=None,
    leaf_fontsize=10,
    figsize='auto',
    gap=0.025,
    show_grid='y',
    title=None,
    out_path=None,
    ax=None,
    **kwargs
):
    '''
    Z : linkage matrix
    labels : list
        List of labels for the leaves of the dendrogram.
    leaf_fontsize : int
        Font size for the labels of the leaves.
    figsize : tuple
        Figure size.
    gap : float
        Proportion of the circle to leave as a "gap" between the dendrogram.
        This gap is placed on the right-hand side of the circle, starting at 0
        degrees (i.e, the horizontal), and puts equal space above and below the
        horizontal.
    show_grid : str
        One of ['x', 'y', True, False].
    title : str
        Title for the plot.
    '''
    def smoothsegment(seg, Nsmooth=100):
        return np.concatenate([[seg[0]], np.linspace(seg[1], seg[2], Nsmooth), [seg[3]]])

    tree = hierarchy.dendrogram(Z, no_plot=True, count_sort=True)

    if kwargs.get('no_plot'):
        pass
    else:
        # 'dcoord' is the width of the branch [???].
        dcoord = np.array(tree['dcoord'])
        dcoord = -np.log(dcoord+1)

        # Rescale icoord: [gap/2, 1-(gap/2)] -> radians; ie, distribute the leaves
        # evenly around the plot.
        # 'icoord' is the leaves and all of the lines parallel to the leaves.
        icoord = np.array(tree['icoord'])
        imax = icoord.max()
        imin = icoord.min()
        # print(f'imin={imin}, imax={imax}')
        icoord = ( (((icoord - imin) / (imax - imin)) * (1-gap)) + gap/2 ) * 2 * np.pi

        if figsize == 'auto':
            figsize = (10, 10)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'})

        # This is the part that makes the actual dendrogram.
        for xs, ys in zip(icoord, dcoord):
            # [TODO] Color the clusters in the dendrogram. Eg, try checking the value
            # of xs,ys to see if it's less than 'color_threshold'. Alternatively,
            # just color the leaves using the 'cluster ids' (don't color the dendrogram at all).
            xs = smoothsegment(xs)
            ys = smoothsegment(ys)
            ax.plot(xs, ys, color="black")

        # Turn off black line around outside of plot.
        ax.spines['polar'].set_visible(False)

        # Put the distance label on the horizontal (0 degrees).
        ax.set_rlabel_position(0)

        if labels:
            n_ticks = len(labels)

            # Set the xtick positions based on the range of icoord, which is in radians.
            imin = icoord.min()
            imax = icoord.max()
            ticks = np.linspace(imin, imax, n_ticks)
            ax.set_xticks(ticks)

            # Match the labels to the tree.
            labels_ = [labels[i] for i in tree['leaves']]
            ax.set_xticklabels(labels_, fontsize=leaf_fontsize)

            # Set the rotation for each label individually.
            gap_in_radians = gap * 2 * np.pi
            start_radians = (gap_in_radians / 2)
            end_radians = (2 * np.pi) - (gap_in_radians / 2)
            radians = np.linspace(start_radians, end_radians, n_ticks)
            radians[np.cos(radians) < 0] = radians[np.cos(radians) < 0] + np.pi
            angles = np.rad2deg(radians)

            # Overwrite the existing plot labels.
            # [TODO] There must be a cleaner way to do this without setting all
            # of the labels first and then re-getting the labels from the figure....
            label_padding = 0.1
            for label, angle in zip(ax.get_xticklabels(), angles):
                x,y = label.get_position()
                lab = ax.text(
                    x,
                    y-label_padding,
                    label.get_text(),
                    transform=label.get_transform(),
                    ha=label.get_ha(),
                    va=label.get_va()
                )
                lab.set_rotation(angle)
            ax.set_xticklabels([])

        # Adjust the grid. The default is to *show* the grid, so we have to
        # explicitly turn it off.
        if not show_grid:
            ax.grid(visible=False)
        elif show_grid == 'y':
            # Show concentric circles. This is the default for this function.
            ax.grid(visible=False, axis='x')
        elif show_grid == 'x':
            ax.grid(visible=False, axis='y')
        else:
            # Show both grids. This is the default in matplotlib.
            pass

        if title is not None:
            ax.set_title(f"Polar= {polar}", fontsize=15)

        savefig(out_path=out_path)

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
        LOGGER.info(f'Saved clusters: {out_path}')

    return clusters


def partition_fullranks(
    path_or_dataframe,
    dendrogram_style='rectangular',
    label_mapper=None,
    out_dendrogram=None,
    out_clusters=None,
    threshold=0,
    **kwargs
):
    '''
    CLI helper.

    Parameters
    ----------
    label_mapper : dict
        Use these labels on the dendrogram. Dict of `{seed: dendrogram_label,
        ... }`.

    Returns
    -------
    dict
    '''
    LOGGER.info('Partitioning....')
    # Z, labels = cluster_hierarchical(path_or_dataframe)
    partition_result = cluster_hierarchical(path_or_dataframe)

    Z = partition_result['linkage']

    labels = partition_result['labels']
    if label_mapper is not None:
        labels = [label_mapper[l] for l in labels]
    elif 'labels' in kwargs:
        labels = kwargs.pop('labels', labels)
    else:
        pass

    if threshold == 'mean':
        threshold = np.mean(Z[:,2])
    elif threshold is None:
        # Default for hierarchy.dendrogram is to use `0.7*max(Z[:,2])` when
        # `color_threshold==None`. Set this explicitely to enable
        # `draw_threshold`, below.
        threshold = 0.7*max(Z[:,2])
    else:
        threshold = float(threshold)

    if dendrogram_style is None:
        # Catch None, bc `None.startswith` raises error.
        tree = {}
    elif dendrogram_style.startswith('r'):
        try:
            tree = plot_dendrogram(
                Z,
                labels=labels,
                color_threshold=threshold,
                out_path=out_dendrogram,
                **kwargs
            )
        except Exception as e:
            LOGGER.error('Plotting failed: %s', str(e))
            tree = {}
    elif dendrogram_style.startswith('p'):
        try:
            tree = plot_dendrogram_polar(
                Z,
                labels=labels,
                out_path=out_dendrogram,
                **kwargs
            )
        except Exception as e:
            LOGGER.error('Plotting failed: %s', str(e))
            tree = {}
    else:
        tree = {}

    clusters = get_clusters(
        Z,
        labels=labels,
        threshold=threshold,
        match_to_leaves=tree.get('leaves'),
        out_path=out_clusters,
    )

    partition_result['tree'] = tree
    partition_result['clusters'] = clusters

    return partition_result


def parse_args(test=None):

    def _valid_threshold_values(arg):
        try:
            if arg == 'mean':
                return arg
            else:
                return float(arg)
        except:
            raise argparse.ArgumentTypeError(f'Threshold must be float or "mean"; you gave "{arg}".')

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
        '--nodetable',
        action='store',
        help='Path to "nodetable" file. This is a TSV file where the first column is the node name (i.e., the seed genes from RWR-fullranks).'
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
        default=0,
        type=_valid_threshold_values,
        help=(
            'Apply threshold to dendrogram. Genes in branches below this threshold will be grouped into clusters; other genes are considered isolates (separate clusters, each with a single gene). Value can be float or "mean". If the value is "mean", then use the mean branch height as the cluster threshold; this can be useful for a first pass.'
        )
    )
    parser.add_argument(
        '--dendrogram-style', '-s',
        action='store',
        choices=['rectangular', 'r', 'polar', 'p', 'none', 'n'],
        default='rectangular',
        help='Plot the dendrogram in rectangular or polar coordinates. Default is rectangular.'
    )
    parser.add_argument(
        '--no-plot',
        action='store_true',
        default=False,
        help='Do not plot the dendrogram.'
    )
    parser.add_argument(
        '--labels-use-names',
        action='store',
        nargs='*',
        type=str,
        help='Label the dendrogram using columns from the nodetable. This is a space-separated list of column names from the nodetable. Pass columns as strings (column names).'
    )
    parser.add_argument(
        '--labels-use-locs',
        action='store',
        nargs='*',
        type=int,
        help='Label the dendrogram using columns from the nodetable. This is a space-separated list of integers indicating columns from the nodetable (0-index, e.g., the first column, which contains the node names, has index 0; the second column has index 1, etc).'
    )
    parser.add_argument(
        '--labels-sep',
        action='store',
        default=' | ',
        help='The separator that will be used if multiple columns from nodetable are used to label the dendrogram.'
    )
    parser.add_argument(
        '--out-dir',
        action='store',
        help='Save dendrogram and clusters to path.'
    )
    parser.add_argument(
        '--out-dendrogram', '-d',
        action='store',
        help='Save dendrogram to path.'
    )
    parser.add_argument(
        '--out-clusters', '-c',
        action='store',
        help='Save clusters to path as tsv file with columns "label", "cluster". When --threshold is 0 (the default) each gene is put into a separate cluster (i.e., every cluster has only a single gene).'
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

    return args


def main():
    args = parse_args()

    # Logging:
    # - Using `LOGGER.setLevel` isn't working, use `logging.basicConfig`
    #   instead.
    # - `logging.basicConfig` should be called *once*.
    # - `logging.basicConfig` also affects settings for imported modules, eg,
    #   matplotlib.
    logger_config = dict(format='[%(asctime)s|%(levelname)s] %(message)s', datefmt='%FT%T')
    if args.verbose == 0:
        # LOGGER.setLevel(logging.WARNING)
        logger_config['level'] = logging.WARNING
    elif args.verbose == 1:
        # LOGGER.setLevel(logging.INFO)
        logger_config['level'] = logging.INFO
    elif args.verbose >= 2:
        # LOGGER.setLevel(logging.DEBUG)
        logger_config['level'] = logging.DEBUG
    logging.basicConfig(**logger_config)

    # LOGGER.debug('debug message')
    # LOGGER.info('info message')
    # LOGGER.warning('warn message')
    # LOGGER.error('error message')
    # LOGGER.critical('critical message')

    if args.out_dir is not None:
        # Use --out-dir with default names, unless another path is explicitely specified.

        if args.out_dendrogram is None:
            # Set the default path for the dendrogram.
            out_dendrogram = os.path.join(args.out_dir, 'dendrogram.png')
        else:
            out_dendrogram = args.out_dendrogram

        if args.out_clusters is None:
            # Set the default path for the clusters.
            out_clusters = os.path.join(args.out_dir, 'clusters.tsv')
        else:
            out_clusters = args.out_clusters
    else:
        out_dendrogram = args.out_dendrogram
        out_clusters = args.out_clusters

    if args.dendrogram_style.startswith('n'):
        dendrogram_style = None
    else:
        dendrogram_style = args.dendrogram_style

    if args.labels_use_names or args.labels_use_locs:
        nodetable = pd.read_table(args.nodetable)
        columns = nodetable.columns.to_list()
        idx_col = columns[0]

        if args.labels_use_locs:
            col_names = [columns[i] for i in args.labels_use_locs]
        else:
            col_names = args.labels_use_names

        label_mapper = nodetable.set_index(idx_col, drop=False)[col_names].agg(args.labels_sep.join, axis=1).to_dict()
    else:
        label_mapper = None


    # Load the fullranks data once.
    fullranks = pd.read_table(args.rwr_fullranks)

    if args.partition:

        result = partition_fullranks(
            path_or_dataframe=fullranks,
            threshold=args.threshold,
            dendrogram_style=dendrogram_style,
            label_mapper=label_mapper,
            out_dendrogram=out_dendrogram,
            out_clusters=out_clusters,
            no_plot=args.no_plot
        )

    return 0


if __name__ == '__main__':
    sys.exit(main())
