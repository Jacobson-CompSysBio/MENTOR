'''
Functional partitioning

This module contains functions for functional partitioning of a network.

[TODO] Refactor:

    1. Make `X` from 'fullranks'.
    2. Apply clustering method to `X`.
    3. Optional: plot dendrogram (only for HC).

References
----------
[1] https://en.wikipedia.org/wiki/Elbow_method_(clustering)
[2] https://en.wikipedia.org/wiki/Root-mean-square_deviation
'''

from functional_partitioning._version import get_version

__version__ = get_version()

import argparse
import os
import sys
import pandas as pd
import numpy as np
import logging
import warnings
import pathlib

from scipy.spatial import distance
from scipy.cluster import hierarchy
from functional_partitioning import _cluster as cluster
from functional_partitioning import _metrics as metrics
from functional_partitioning import _rwrtoolkit as rwrtoolkit
from functional_partitioning import _plot as plot


LOGGER = logging.getLogger(__name__)


def parse_args(test=None):

    def _valid_threshold_values(arg):
        try:
            if isinstance(arg, str):
                arg = arg.lower()

            if arg is None:
                return None
            elif arg == 'mean':
                return arg
            elif arg == 'best_chi':
                return arg
            elif arg == 'none':
                return None
            else:
                return float(arg)
        except:
            raise argparse.ArgumentTypeError(f'Threshold must be one of "mean", "best_chi", or <float>; you gave "{arg}".')

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
        help='Perform functional partitioning on "seed genes" from RWR fullranks file. This is the default.'
    )
    parser.add_argument(
        '--no-partition',
        action='store_false',
        dest='partition',
        help='Do not perform functional partitioning.'
    )
    parser.add_argument(
        '--threshold', '-t',
        action='store',
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
        help='Plot the dendrogram in rectangular or polar coordinates. If "none", then do not plot the dendrogram (this is redundant with --no-plot).'
    )
    parser.add_argument(
        '--no-plot',
        action='store_true',
        default=False,
        help='Do not plot the dendrogram.'
    )
    parser.add_argument(
        '--labels-use-clusters',
        action='store_true',
        default=False,
        help=''
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
        '--outdir',
        action='store',
        type=pathlib.Path,
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
        '--path-to-conda-env',
        action='store',
        help=''
    )
    parser.add_argument(
        '--path-to-rwrtoolkit',
        action='store',
        help=''
    )
    parser.add_argument(
        '--multiplex',
        action='store',
        help=''
    )
    parser.add_argument(
        '--geneset',
        action='store',
        help=''
    )
    parser.add_argument(
        '--method',
        action='store',
        default='singletons',
        help=''
    )
    parser.add_argument(
        '--folds',
        action='store',
        help=''
    )
    parser.add_argument(
        '--restart',
        action='store',
        help=''
    )
    parser.add_argument(
        '--tau',
        action='store',
        help=''
    )
    parser.add_argument(
        '--numranked',
        action='store',
        help=''
    )
    parser.add_argument(
        '--modname',
        action='store',
        help=''
    )
    parser.add_argument(
        '--plot',
        action='store',
        help=''
    )
    parser.add_argument(
        '--threads',
        action='store',
        help=''
    )
    parser.add_argument(
        '--init-test-fullranks',
        action='store',
        help='Create fullranks file for testing at the given path.'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='count',
        default=0,
        help='Default: WARNING; once: INFO; twice: DEBUG'
    )
    parser.add_argument(
        '--version',
        action='store_true',
        default=False,
        help='Print version and exit.'
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
    # print(args)

    # Logging:
    # - `LOGGER.setLevel` isn't working, use `logging.basicConfig` instead.
    # - `logging.basicConfig` should be called *once*.
    # - `logging.basicConfig` also affects settings for imported modules, eg, matplotlib.
    logger_config = dict(
        format='[%(asctime)s|%(levelname)s] %(message)s',
        datefmt='%FT%T',
        level=logging.WARNING
    )
    if args.verbose == 1:
        logger_config['level'] = logging.INFO
    elif args.verbose >= 2:
        logger_config['level'] = logging.DEBUG
    logging.basicConfig(**logger_config)

    # Test logging messages.
    # LOGGER.debug('debug message')
    # LOGGER.info('info message')
    # LOGGER.warning('warn message')
    # LOGGER.error('error message')
    # LOGGER.critical('critical message')

    if args.version:
        print(__version__)
        sys.exit(0)

    # Create dummy data for testing.
    if args.init_test_fullranks:
        from functional_partitioning import _datasets as datasets
        fullranks = datasets.make_fullranks_table()
        # print('fullranks:')
        # print(fullranks)
        fullranks.to_csv(args.init_test_fullranks, sep='\t', index=False)

    # Use --out-dir with default names, unless another path is explicitely
    # specified.
    if args.out_clusters is not None:
        # Set the default path for the clusters.
        out_clusters = args.out_clusters
    elif args.outdir is not None:
        out_clusters = os.path.join(args.outdir, 'clusters.tsv')
    else:
        out_clusters = None

    # Use --out-dir with default names, unless another path is explicitely
    # specified.
    if args.no_plot:
        out_dendrogram = None
    elif args.out_dendrogram is not None:
        # Set the default path for the dendrogram.
        out_dendrogram = args.out_dendrogram
    elif args.outdir is not None:
        out_dendrogram = os.path.join(args.outdir, 'dendrogram.png')
    else:
        out_dendrogram = None

    # Set dendrogram style (rectangular or polar).
    if out_dendrogram is not None and args.dendrogram_style.startswith(('r', 'p')):
        dendrogram_style = args.dendrogram_style
    else:
        dendrogram_style = None

    # Run RWR-singletons or get the fullranks file.
    if args.multiplex and args.geneset:
        # Run RWR-singletons.
        command = rwrtoolkit.rwr_singletons(
            path_to_conda_env=args.path_to_conda_env,
            path_to_rwrtoolkit=args.path_to_rwrtoolkit,
            data=args.multiplex,
            geneset=args.geneset,
            method=args.method,
            folds=args.folds,
            restart=args.restart,
            tau=args.tau,
            numranked=args.numranked,
            outdir=args.outdir,
            modname=args.modname,
            plot=args.plot,
            threads=args.threads,
            verbose=args.verbose
        )
        # print(command)
        res = rwrtoolkit.run(command)
        if res['returncode'] != 0:
            LOGGER.error('RWR-singletons failed.')
            LOGGER.error(command)
            LOGGER.error(res.stderr)
            sys.exit(1)
        # print(res)

        rwrtoolkit.compress_results(args.outdir)

        try:
            path_to_fullranks = next(args.outdir.glob('RWR*fullranks*'))
            # print(path_to_fullranks)
        except StopIteration:
            LOGGER.error('Cannot find fullranks file.')
            sys.exit(1)

    else:
        # Read the fullranks file from existing RWR-singletons results.
        path_to_fullranks = args.rwr_fullranks

    # Run functional partitioning or exit.
    if args.partition:
        # Run functional partitioning.
        X_ranks, X_scores, labels, max_rank = rwrtoolkit.fullranks_to_matrix(
            path_to_fullranks,
            max_rank='elbow',
            drop_missing=True
        )

        linkage_matrix = cluster.cluster_hierarchical(
            X_ranks.fillna(0),
            corr_method='spearman',
            linkage_method='average'
        )
        threshold = cluster.calc_threshold(
            linkage_matrix,
            args.threshold,
            scores=X_scores.fillna(X_scores.max(axis=1))
        )
        # print('threshold:', threshold)

        # print('linkage matrix:')
        # print(linkage_matrix)
        clusters = cluster.get_clusters(
            linkage_matrix,
            labels=labels,
            threshold=threshold,
            n_clusters=None,
            match_to_leaves=None,
            out_path=out_clusters
        )
        # print('clusters:')
        # print(clusters)

        if args.labels_use_clusters:
            label_mapper = plot.make_label_mapper(
                nodetable=clusters,
                use_locs=[0, 1], # List.
                sep=args.labels_sep
            )
            labels = [label_mapper.get(l, l) for l in labels]
        elif args.labels_use_names or args.labels_use_locs:
            label_mapper = plot.make_label_mapper(
                nodetable=args.nodetable,
                use_names=args.labels_use_names,
                use_locs=args.labels_use_locs,
                sep=args.labels_sep
            )
            labels = [label_mapper.get(l, l) for l in labels]
        else:
            label_mapper = None
        # print(label_mapper)
        # print(labels)

        if dendrogram_style is None:
            # Catch None, bc `None.startswith` raises error.
            tree = {}
        elif dendrogram_style.startswith('r'):
            # Rectangular dendrogram.
            try:
                tree = plot.plot_dendrogram(
                    linkage_matrix,
                    labels=labels,
                    color_threshold=threshold,
                    out_path=out_dendrogram,
                    no_plot=args.no_plot
                )
            except Exception as e:
                # LOGGER.error('Plotting failed: %s', str(e))
                warnings.warn('[WARNING] Unable to draw the dendrogram, see error message:\n %s' % str(e))
                tree = {}
        elif dendrogram_style.startswith('p'):
            # Polar dendrogram.
            try:
                tree = plot.plot_dendrogram_polar(
                    linkage_matrix,
                    labels=labels,
                    out_path=out_dendrogram,
                    no_plot=args.no_plot
                )
            except Exception as e:
                # LOGGER.error('Plotting failed: %s', str(e))
                warnings.warn('[WARNING] Unable to draw the dendrogram, see error message:\n %s' % str(e))
                tree = {}
        else:
            tree = {}
    else:
        # Exit.
        pass

    return 0

# END.
