
from mentor._version import __version__
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
from mentor import _cluster as cluster
from mentor import _metrics as metrics
from mentor import _rwrtoolkit as rwrtoolkit
from mentor import _fancydend as fd

LOGGER = logging.getLogger(__name__)

def parse_args(test=None):

    parser = argparse.ArgumentParser(
        description='Partition seeds from `RWR-CV --method=singletons ...` into clusters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--outdir',
        action='store',
        type=pathlib.Path,
        help='Location to save dendrogram and clusters'
    )
    parser.add_argument(
        '--outfile',
        action='store',
        help='Description to append to filenames'
    )
    parser.add_argument(
        '--multiplex',
        action='store',
        help='Full path to multiplex network'
    )
    parser.add_argument(
        '--geneset',
        action='store',
        help='Full path to geneset'
    )
    parser.add_argument(
        '--threads',
        action='store',
        help=''
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
        help='Print version and exit'
    )
    parser.add_argument(
        '--distances',
        action='store',
        help='Full path to the distance matrix file'
    )
    parser.add_argument(
        '--clusters',
        action='store',
        default=3,
        type=int,
        help='Number of initial clusters desired in dendrogram (inital # if using --subcluster)'
    )
    parser.add_argument(
        '--map',
        action='store',
        help='Full path to the gene symbol mapping file'
    )
    parser.add_argument(
        '--subcluster',
        action='store_true',
        default=False,
        help='Specify if you want to subcluster dendrogram'
    )
    parser.add_argument(
        '--increment',
        action='store',
        default = 5,
        type=int,
        help='Increment to use in subclustering'
    )
    parser.add_argument(
        '--maxsize',
        action='store',
        default = 40,
        type=int,
        help='Maximum size for clades if subclustering'
    )
    parser.add_argument(
        '--heatmaps',
        action='store',
        help='Full path to the heatmap file if you want to include a heatmap'
    )
    parser.add_argument(
        '--plottype',
        action='store',
        default='rectangular',
        help='Type of dendrogram (either rectangular or polar)'
    )
    parser.add_argument(
        '--reordercols',
        action='store_true',
        default=False,
        help='Specify if you want to perform hierarchical clustering on the columns of the heatmap'
    )
    parser.add_argument(
        '--legendtitle',
        action='store',
        default="Value,Group",
        help='Titles to give to legends; for rectangular dendrogram title to give to continuous legend; for polar dendrogram titles to give to continuous legend and factor legend'
    )
    parser.add_argument(
        '--squish',
        action='store',
        help='If continuous columns in heatmap present squish the upper and lower bounds of color scheme to lower_bound,upper_bound'
    )
    parser.add_argument(
        '--relwidths',
        action='store',
        default='1,1',
        help='If including a heatmap then adjust the relative widths of the dendrogram to the heatmap with dend_width,heatmap_width; only used for rectangular dendrogram'
    )
    parser.add_argument(
        '--clusterlabelsize',
        action='store',
        default=2.5,
        help='Size of cluster labels associated with polar dendrogram'
    )
    parser.add_argument(
        '--heatmaplabelsize',
        action='store',
        default=0.75,
        help='Size of the labels associated with the polar dendrogram heatmap (gene names)'
    )
    parser.add_argument(
        '--groupcolors',
        action='store',
        default=None,
        help='Colors for the group blocks in polar dendrogram'
    )
    parser.add_argument(
        '--trackheight',
        action='store',
        default='0.2,0.2,0.2',
        help='Width of the tracks in the polar dendrogram (heatmap1,heatmap2,dendrogram)'
    )
    parser.add_argument(
        '--highlightindex',
        action='store',
        default=None,
        help='Sector(s) to highlight on polar dendrogram (sector1,sector2,...,sectorN)'
    )
    parser.add_argument(
        '--highlightcolor',
        action='store',
        default='#34EBDC',
        help='Color to use for highlighted sectors on polar dendrogram'
    )
    parser.add_argument(
        '--plotwidth',
        action='store',
        default=30,
        type=int,
        help='Width of the plot'
    )
    parser.add_argument(
        '--plotheight',
        action='store',
        default=None,
        type=int,
        help='Height of the plot'
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
    logger_config = dict(
        format = '[%(asctime)s|%(levelname)s] %(message)s',
        datefmt = '%FT%T',
        level = logging.WARNING
    )
    if args.verbose == 1:
        logger_config['level'] = logging.INFO
    elif args.verbose >= 2:
        logger_config['level'] = logging.DEBUG
    logging.basicConfig(**logger_config)
    if args.version:
        print(__version__)
        sys.exit(0)
    LOGGER.debug(args)
    if args.outfile is not None:
        diss_mat = args.outfile + '-dissimilarity-matrix.tsv'
        diss_stats = args.outfile + '-dissimilarity-stats.tsv'
        diss_dist = args.outfile + '-distribution-of-pairwise-dissimilarities.png'
    else:
        diss_mat = 'dissimilarity-matrix.tsv'
        diss_stats = 'dissimilarity-stats.tsv'
        diss_dist = 'distribution-of-pairwise-dissimilarities.png'
    if args.outdir is not None:
        out_dissimilarity_matrix = args.outdir / diss_mat
        out_dissimilarity_stats = args.outdir / diss_stats
        out_dissimilarity_distribution = args.outdir / diss_dist
    else:
        out_dissimilarity_matrix = None
        out_dissimilarity_stats = None
        out_dissimilarity_distribution = None
    if args.multiplex and args.geneset:
        command = rwrtoolkit.rwr_singletons(
            data = args.multiplex,
            geneset = args.geneset,
            outdir = args.outdir,
            threads = args.threads,
            verbose = args.verbose
        )
        print("\nrunning RWRtoolkit singletons command: " + command)
        res = rwrtoolkit.run(command)
        if res['returncode'] != 0:
            LOGGER.error('RWR-singletons failed.')
            LOGGER.error(command)
            sys.exit(1)
        rwrtoolkit.compress_results(args.outdir)
        try:
            path_to_fullranks = next(args.outdir.glob('RWR*fullranks*'))
        except StopIteration:
            LOGGER.error('Cannot find fullranks file.')
            sys.exit(1)
        print("\nrunning MENTOR")
        X, labels = rwrtoolkit.transform_fullranks(
            path_to_fullranks,
            drop_missing = True,
            max_rank = 'elbow',
        )
        print("\nrunning clustering")
        mod = cluster.HierarchicalClustering(
            metric = metrics.spearman_d,
            memory = None,
            connectivity = None,
            linkage = "average",
            compute_distances = False,
            compute_linkage_matrix = True
        )
        mod.fit(X)
        if out_dissimilarity_matrix is not None:
            print("\nsaving dissimilarity matrix")
            out_dissimilarity_matrix.parent.mkdir(parents = False,exist_ok = True)
            dmat = pd.DataFrame(
                distance.squareform(mod.pairwise_distances,checks = False),
                index = labels,
                columns = labels
            )
            dmat.to_csv(out_dissimilarity_matrix,sep = '\t')
            fd.fancy_dendrogram(
                distances = out_dissimilarity_matrix,
                clusters = args.clusters,
                map = args.map,
                outdir = args.outdir,
                outfile = args.outfile,
                subcluster = args.subcluster,
                increment = args.increment,
                maxsize = args.maxsize,
                heatmaps = args.heatmaps,
                plottype = args.plottype,
                reordercols = args.reordercols,
                legendtitle = args.legendtitle,
                squish = args.squish,
                relwidths = args.relwidths,
                clusterlabelsize = args.clusterlabelsize,
                heatmaplabelsize = args.heatmaplabelsize,
                groupcolors = args.groupcolors,
                trackheight = args.trackheight,
                highlightindex = args.highlightindex,
                highlightcolor = args.highlightcolor,
                plotwidth = args.plotwidth,
                plotheight = args.plotheight
            )
        else:
            dmat = None
        if out_dissimilarity_stats is not None:
            print("\nsaving dissimilarity matrix statistics")
            out_dissimilarity_stats.parent.mkdir(parents = False,exist_ok = True)
            with open(out_dissimilarity_stats,'w') as f:
                try:
                    dist_summary = metrics.summarize_pairwise_dissimilarities(mod.pairwise_distances,mod.labels_)
                    for key,value in dist_summary.items():
                        print(key,value,sep = '\t',file = f)
                except:
                    pass
                try:
                    chi = metrics.calinski_harabasz_score_(X,mod.labels_)
                    print('calinski_harabasz_score',chi,sep = '\t',file = f)
                except:
                    pass
    elif args.distances is not None:
        fd.fancy_dendrogram(
            distances = args.distances,
            clusters = args.clusters,
            map = args.map,
            outdir = args.outdir,
            outfile = args.outfile,
            subcluster = args.subcluster,
            increment = args.increment,
            maxsize = args.maxsize,
            heatmaps = args.heatmaps,
            plottype = args.plottype,
            reordercols = args.reordercols,
            legendtitle = args.legendtitle,
            squish = args.squish,
            relwidths = args.relwidths,
            clusterlabelsize = args.clusterlabelsize,
            heatmaplabelsize = args.heatmaplabelsize,
            groupcolors = args.groupcolors,
            trackheight = args.trackheight,
            highlightindex = args.highlightindex,
            highlightcolor = args.highlightcolor,
            plotwidth = args.plotwidth,
            plotheight = args.plotheight
        )
    else:
        pass
    return 0
