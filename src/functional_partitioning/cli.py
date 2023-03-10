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

from functional_partitioning._version import __version__

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
from functional_partitioning import _utils as utils


LOGGER = logging.getLogger(__name__)


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
        '--cut-threshold', '-t',
        action='store',
        default=0.3,
        type=float,
        help=('Cut the dendrogram at this threshold. Only used if `--cut-method=hard`.')
    )
    parser.add_argument(
        '--cut-method', '-m',
        action='store',
        choices=['dynamic', 'hard', 'none'],
        default='dynamic',
        help=(f'If `dynamic`, use dynamicTreeCut to determine clusters without a hard threshold. If `none`, return all the clusterings. Otherwise, use with `--cut-threshold` to provide a numeric cut threshold.')
    )
    parser.add_argument(
        '--dendrogram-style', '-s',
        action='store',
        choices=['rectangular', 'r', 'polar', 'p', 'none', 'n'],
        default='rectangular',
        help='Plot the dendrogram in rectangular or polar coordinates. If "none", then do not plot the dendrogram (this is redundant with --no-plot).'
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
        '--no-plot',
        action='store_true',
        default=False,
        help='Do not plot the dendrogram.'
    )
    parser.add_argument(
        '--out-clusters', '-c',
        action='store',
        help='Save clusters to path as tsv file with columns "label", "cluster". When --threshold is 0 (the default) each gene is put into a separate cluster (i.e., every cluster has only a single gene).'
    )
    parser.add_argument(
        '--no-clusters',
        action='store_true',
        default=False,
        help='Do not export clusters to file.'
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

    LOGGER.debug(args)

    # Create dummy data for testing.
    if args.init_test_fullranks:
        from functional_partitioning import _datasets as datasets
        fullranks = datasets.make_fullranks_table()
        # print('fullranks:')
        # print(fullranks)
        fullranks.to_csv(args.init_test_fullranks, sep='\t', index=False)

    # Use --out-dir with default names, unless another path is explicitely
    # specified.
    if args.no_clusters:
        out_clusters = None
        out_dissimilarity_matrix = None
        out_dissimilarity_stats = None
    elif args.out_clusters is not None:
        out_clusters = args.out_clusters
        out_dissimilarity_matrix = args.out_clusters.with_name('dissimilarity-matrix.tsv')
        out_dissimilarity_stats = args.out_clusters.with_name('dissimilarity-stats.tsv')
    elif args.outdir is not None:
        # Set the default path for the clusters.
        out_clusters = args.outdir / 'clusters.tsv'
        out_dissimilarity_matrix = args.outdir / 'dissimilarity-matrix.tsv'
        out_dissimilarity_stats = args.outdir / 'dissimilarity-stats.tsv'
    else:
        out_clusters = None
        out_dissimilarity_matrix = None
        out_dissimilarity_stats = None

    # Use --out-dir with default names, unless another path is explicitely
    # specified.
    if args.no_plot:
        out_dendrogram = None
        out_dissimilarity_distribution = None
    elif args.out_dendrogram is not None:
        out_dendrogram = args.out_dendrogram
        out_dissimilarity_distribution = args.out_dendrogram.with_name('distribution-of-pairwise-dissimilarities.png')
    elif args.outdir is not None:
        # Set the default path for the dendrogram.
        out_dendrogram = args.outdir / 'dendrogram.png'
        out_dissimilarity_distribution = args.outdir / 'distribution-of-pairwise-dissimilarities.png'
    else:
        out_dendrogram = None
        out_dissimilarity_distribution = None

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

        res = rwrtoolkit.run(command)
        if res['returncode'] != 0:
            LOGGER.error('RWR-singletons failed.')
            LOGGER.error(command)
            LOGGER.error(res.stderr)
            sys.exit(1)

        rwrtoolkit.compress_results(args.outdir)

        try:
            path_to_fullranks = next(args.outdir.glob('RWR*fullranks*'))
        except StopIteration:
            LOGGER.error('Cannot find fullranks file.')
            sys.exit(1)

    else:
        # Read the fullranks file from existing RWR-singletons results.
        path_to_fullranks = args.rwr_fullranks

    # Run functional partitioning or exit.
    if args.partition:
        # Run functional partitioning.
        scores, ranks, labels = rwrtoolkit.fullranks_to_matrix(
            path_to_fullranks,
            drop_missing=True,
            verify=True
        )

        X = utils.get_top_ranked(scores, ranks=ranks, max_rank='elbow')

        if args.cut_method == 'dynamic':
            cut_method = 'cutreeHybrid'
            cut_threshold = args.cut_threshold
        elif args.cut_method == 'none':
            cut_method = None
            cut_threshold = None
        else:
            cut_method = args.cut_method
            cut_threshold = args.cut_threshold

        mod = cluster.HierarchicalClustering(
            n_clusters=None,
            metric=metrics.spearman_d,
            cut_method=cut_method,
            cut_threshold=cut_threshold,
            memory=None,
            connectivity=None,
            compute_full_tree="auto",
            linkage="average",
            compute_distances=False,
            compute_linkage_matrix=True,
            compute_dendrogram=True,
        )
        mod.fit(X)
        clusters = pd.DataFrame(mod.labels_, index=labels)
        threshold = mod.cut_threshold_

        if out_clusters is not None:
            clusters.to_csv(out_clusters, sep='\t')
            LOGGER.info(f'Clusters saved to {out_clusters}')

        if out_dissimilarity_matrix is not None:
            dmat = pd.DataFrame(
                distance.squareform(mod.pairwise_distances, checks=False),
                index=labels,
                columns=labels
            )
            dmat.to_csv(out_dissimilarity_matrix, sep='\t')
            LOGGER.info(f'dissimilarity matrix saved to {out_dissimilarity_matrix}')
        else:
            dmat = None

        if out_dissimilarity_stats is not None:
            with open(out_dissimilarity_stats, 'w') as f:
                # print(
                #     'mean of pair-wise dissimilarities',
                #     np.mean(mod.pairwise_distances),
                #     sep='\t',
                #     file=f
                # )
                # print(
                #     'median of pair-wise dissimilarities',
                #     np.median(mod.pairwise_distances),
                #     sep='\t',
                #     file=f
                # )
                try:
                    dist_summary = metrics.summarize_pairwise_dissimilarities(
                        mod.pairwise_distances,
                        mod.labels_
                    )
                    for key, value in dist_summary.items():
                        print(
                            key,
                            value,
                            sep='\t',
                            file=f
                        )
                except:
                    pass
                try:
                    chi = metrics.calinski_harabasz_score_(
                        X,
                        mod.labels_
                    )
                    print(
                        'calinski_harabasz_score',
                        chi,
                        sep='\t',
                        file=f
                    )
                except:
                    pass


        if out_dissimilarity_distribution is not None:
            # print(mod.pairwise_distances)
            # print(np.mean(mod.pairwise_distances))
            plot.pairwise_distances_violin(
                mod.pairwise_distances,
                out_path=out_dissimilarity_distribution
            )

        if out_dendrogram is not None:
            # Set up the leaf labels for the dendrogram.
            if args.labels_use_clusters:
                # Label the leaves with the cluster ID.
                label_mapper = plot.make_label_mapper(
                    nodetable=clusters,
                    use_locs=[0, 1], # List.
                    sep=args.labels_sep
                )
                labels = [label_mapper.get(l, l) for l in labels]
            elif args.labels_use_names or args.labels_use_locs:
                # Label the leaves using the node table.
                label_mapper = plot.make_label_mapper(
                    nodetable=args.nodetable,
                    use_names=args.labels_use_names,
                    use_locs=args.labels_use_locs,
                    sep=args.labels_sep
                )
                labels = [label_mapper.get(l, l) for l in labels]
            else:
                # Use default labels.
                pass

            tree = plot.draw_dendrogram(
                dendrogram_style=dendrogram_style,
                linkage_matrix=mod.linkage_matrix,
                labels=labels,
                threshold=threshold,
                out_path=out_dendrogram,
                no_plot=args.no_plot
            )
    else:
        # Exit.
        pass

    return 0


# END.
