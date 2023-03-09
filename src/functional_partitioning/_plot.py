'''
Plotting functions.
'''

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
        # The default for `read_table` is to set a numeric index, ie, all of
        # the data will appear in columns and will not be used as the index.
        nodetable = pd.read_table(nodetable, index_col=0)
    else:
        nodetable = nodetable.copy()

    # Move the index to a column, preserving the 'name' of the index.
    # nodetable.insert(0, '__index__', nodetable.index)
    idx = nodetable.index
    nodetable = nodetable.reset_index()
    nodetable.index = idx

    columns = nodetable.columns.to_list()
    if use_locs:
        col_names = [columns[i] for i in use_locs]
    elif use_cols:
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


def plot_dendrogram(
    Z,
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
    if kwargs.get('color_threshold') is None:
        # Do not allow 'color_threshold' to be None.
        # - Set it to 0 to prevent dendrogram from using default color_threshold.
        # - Protect against error: '>' not supported between instances of 'NoneType' and 'int'
        kwargs['color_threshold'] = 0

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
        if title is not None:
            plt.gca().set_title(f"{title}", fontsize=15)
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
            ax.set_title(f"{title}", fontsize=15)

        savefig(out_path=out_path)

    return tree


# END.
