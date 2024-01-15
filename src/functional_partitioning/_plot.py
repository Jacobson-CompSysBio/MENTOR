'''
Plotting functions.
'''

import logging
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings

from scipy.cluster import hierarchy

DPI = 300

LOGGER = logging.getLogger(__name__)

def make_label_mapper(nodetable=None, use_names=False, use_locs=False, sep=' | '):
    '''
    Create a label mapper.

    Parameters
    ----------
    nodetable : {pandas.DataFrame, str}
        A pandas DataFrame or path to a tsv file. The index should match the
        node labels from the input data.
    use_names : list
        List of column names in `nodetable` to use as labels.
    use_locs : list
        List of ints corresponding to column locations in `nodetable` to use as
        labels.
    sep : str
        Separator to use when joining multiple columns.

    Returns
    -------
    dict
        {<node name> : <label> , ...}
    '''
    if not use_names and not use_locs:
        raise ValueError('You must provide at least one of `use_names` or `use_locs`.')

    def join(x, sep=sep):
        return sep.join([str(i) for i in x])

    if isinstance(nodetable, str):
        nodetable = pd.read_table(nodetable, index_col=0)
    else:
        nodetable = nodetable.copy()

    idx = nodetable.index
    nodetable = nodetable.reset_index(col_level=0)
    nodetable.index = idx

    columns = nodetable.columns.to_list()
    if use_locs:
        col_names = [columns[i] for i in use_locs]
    elif use_names:
        col_names = use_names
    label_mapper = nodetable[col_names].agg(join, axis=1).to_dict()

    return label_mapper


def savefig(out_path=None, **kwargs):
    kwargs.setdefault('dpi', DPI)
    kwargs.setdefault('bbox_inches', 'tight')
    if out_path is not None:
        try:
            plt.savefig(out_path, **kwargs)
            LOGGER.info('Saved figure: %s', out_path)
        except Exception as e:
            LOGGER.error('Failed to save figure: %s', str(e))


def _plot_dendrogram_rectangular(
    linkage_matrix=None,
    out_path=None,
    figsize='auto',
    draw_threshold=True,
    title=None,
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
    _orientation = kwargs.get('orientation', 'left')

    if kwargs.get('color_threshold') is None:
        kwargs['color_threshold'] = 0

    if kwargs.get('no_plot'):
        pass
    elif plt.get_fignums() and kwargs.get('ax') is None:
        kwargs['ax'] = plt.gca()
    elif kwargs.get('ax') is None:
        if figsize == 'auto':
            width = 5
            height = np.shape(linkage_matrix)[0] * 0.2
            if height < 10:
                height = 10
            if _orientation == 'top' or _orientation == 'bottom':
                width, height = height, width
            figsize_ = (width, height)
        else:
            figsize_ = figsize
        plt.rc('figure', facecolor='white')
        plt.figure(figsize=figsize_)

    if kwargs.get('no_plot'):
        pass
    elif draw_threshold and kwargs.get('color_threshold', 0) > 0:
        if _orientation == 'left' or _orientation == 'right':
            plot_line = plt.axvline
        elif _orientation == 'top' or _orientation == 'bottom':
            plot_line = plt.axhline
            if figsize == 'auto':
                width, height = height, width
                figsize_ = (width, height)
        else:
            raise ValueError(f'`orientation` must be one of ["top", "bottom", "left", "right"]: {_orientation}')
        plot_line(kwargs.get('color_threshold'), c='k', linewidth=1, linestyle='dotted')

    tree = hierarchy.dendrogram(
        linkage_matrix,
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
        if title is not None:
            plt.gca().set_title(f"{title}", fontsize=15)
        if out_path:
            savefig(out_path=out_path)
            plt.close()

    return tree


def _plot_dendrogram_polar(
    linkage_matrix=None,
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
    linkage_matrix : numpy.ndarray
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

    tree = hierarchy.dendrogram(linkage_matrix, no_plot=True, count_sort=True)

    if kwargs.get('no_plot'):
        pass
    else:
        dcoord = np.array(tree['dcoord'])
        dcoord = -np.log(dcoord+1)
        icoord = np.array(tree['icoord'])
        imax = icoord.max()
        imin = icoord.min()
        icoord = ( (((icoord - imin) / (imax - imin)) * (1-gap)) + gap/2 ) * 2 * np.pi

        if figsize == 'auto':
            figsize = (10, 10)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'})

        # this is the part that makes the actual dendrogram
        for xs, ys in zip(icoord, dcoord):
            xs = smoothsegment(xs)
            ys = smoothsegment(ys)
            ax.plot(xs, ys, color="black")

        ax.spines['polar'].set_visible(False)

        ax.set_rlabel_position(0)

        if labels:
            n_ticks = len(labels)
            imin = icoord.min()
            imax = icoord.max()
            ticks = np.linspace(imin, imax, n_ticks)
            ax.set_xticks(ticks)
            labels_ = [labels[i] for i in tree['leaves']]
            ax.set_xticklabels(labels_, fontsize=leaf_fontsize)
            gap_in_radians = gap * 2 * np.pi
            start_radians = (gap_in_radians / 2)
            end_radians = (2 * np.pi) - (gap_in_radians / 2)
            radians = np.linspace(start_radians, end_radians, n_ticks)
            radians[np.cos(radians) < 0] = radians[np.cos(radians) < 0] + np.pi
            angles = np.rad2deg(radians)
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

        if not show_grid:
            ax.grid(visible=False)
        elif show_grid == 'y':
            ax.grid(visible=False, axis='x')
        elif show_grid == 'x':
            ax.grid(visible=False, axis='y')
        else:
            pass

        if title is not None:
            ax.set_title(f"{title}", fontsize=15)

        if out_path:
            savefig(out_path=out_path)
            plt.close()

    return tree


def draw_dendrogram(dendrogram_style=None, **kwargs):
    tree = {}
    if dendrogram_style and dendrogram_style.startswith('r'):
        LOGGER.info('Drawing rectangular dendrogram.')
        try:
            tree = _plot_dendrogram_rectangular(
                linkage_matrix=kwargs.get('linkage_matrix'),
                labels=kwargs.get('labels'),
                color_threshold=kwargs.get('threshold'),
                out_path=kwargs.get('out_path'),
                no_plot=kwargs.get('no_plot')
            )
        except Exception as e:
            warnings.warn(
                'Unable to draw the dendrogram, see error message:\n %s' % str(e),
                UserWarning
            )
            tree = {}
    elif dendrogram_style and dendrogram_style.startswith('p'):
        LOGGER.info('Drawing polar dendrogram.')
        try:
            tree = _plot_dendrogram_polar(
                linkage_matrix=kwargs.get('linkage_matrix'),
                labels=kwargs.get('labels'),
                out_path=kwargs.get('out_path'),
                no_plot=kwargs.get('no_plot')
            )
        except Exception as e:
            warnings.warn(
                'Unable to draw the dendrogram, see error message:\n %s' % str(e),
                UserWarning
            )
            tree = {}
    else:
        LOGGER.debug('No dendrogram requested.')

    return tree

def pairwise_distances_violin(
    data,
    title=None,
    out_path=None,
    ax=None,
    ):

    if ax is None:
        plt.figure()
        ax = plt.gca()
    p = sns.violinplot(
        data=data,
        orient='h',
        cut=0,
        ax=ax
    )
    plt.xlabel('dissimilarity')
    plt.ylabel('distribution')
    if title:
        plt.title(title)
    if out_path:
        savefig(out_path=out_path)
        plt.close()
