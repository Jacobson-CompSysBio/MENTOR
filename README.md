% Functional Partitioning
% J. Izaak Miller

<!-- BADGES {{{ -->
[![DOI](https://zenodo.org/badge/526587672.svg)](https://zenodo.org/badge/latestdoi/526587672)
<!-- BADGES }}} -->

MENTOR: Mechanistic Exploration of Networks for Team-based Omics Research
=======================

- Integrate multi-omics data using multiplex networks.
- Identify functionally-related groups of genes using random walk with restart
  on multiplex networks.

Installation
============

Standard installation:

```sh
$ git clone https://github.com/izaakm/jail-functional-partitioning
$ cd jail-functional-partitioning
$ make
```

If you prefer to install dependencies with conda, you can use the provided
`environment.yml` file:

```sh
$ git clone https://github.com/izaakm/jail-functional-partitioning
$ cd jail-functional-partitioning
$ conda env update -f ./environment.yml
$ make
```

Or create a new environment for 'functional-partitioning':

```sh
$ git clone https://github.com/izaakm/jail-functional-partitioning
$ cd jail-functional-partitioning
$ conda env create -f ./environment.yml
$ conda activate functional-partitioning
$ make
```

Verify the installation:

```sh
$ functional_partitioning --help
```


Usage
=====

```

usage: mentor [-h] [--rwr-fullranks RWR_FULLRANKS] [--partition] [--outdir OUTDIR] [--path-to-conda-env CONDA_ENV] [--path-to-rwrtoolkit RWRTOOLKIT] [--multiplex MULTIPLEX] [--geneset GENESET] [--method METHOD] [--folds FOLDS] [--restart RESTART] [--tau TAU] [--numranked NUMRANKED] [--modname MODNAME] [--threads THREADS] [--verbose] [--version] [--distances DISTANCES] [--clusters CLUSTERS] [--map MAP] [--subcluster] [--increment INCREMENT] [--maxsize MAXSIZE] [--heatmaps HEATMAPS] [--pcutoff PCUTOFF] [--squish LOWER,UPPER] [--relwidths DEND,HEAT] [--plotwidth PLOTWIDTH]

arguments:
  -h, --help                           Show this help message and exit.
  --rwr-fullranks RWR_FULLRANKS        Path to "fullranks" file from `RWR-CV--method=singletons`
                                       (default: None).
  --partition, -p                      Perform functional partitioning on "seed genes" from RWR
                                       fullranks file (default: False).
  --outdir OUTDIR                      Save output to path (default: None).
  --path-to-conda-env CONDA_ENV        Path to conda environment.
  --path-to-rwrtoolkit RWRTOOLKIT      Path to RWRToolkit.
  --multiplex MULTIPLEX                Path to multiplex network.
  --geneset GENESET                    Path to gene set file. 
  --method METHOD                      Method for RWR-CV (default: 'singletons').
  --folds FOLDS                        Folds for RWR-CV.
  --restart RESTART                    Restart for RWR-CV (default: 0.7). 
  --tau TAU                            Tau for RWR-CV.
  --numranked NUMRANKED                
  --modname MODNAME
  --threads THREADS
  --verbose, -v                        Default: WARNING; once: INFO; twice: DEBUG (default: 0)
  --version                            Print version and exit (default: False).
  --distances DISTANCES                Path to dissimilarity-matrix.tsv (default: None).
  --clusters CLLUSTERS                 Number of clusters for dendrogram (default: 10).
  --map                                Path to gene ensembl ID mapping file (default: None).
  --subcluster                         Subcluster the dendrogram (default: False).   
  --increment INCREMENT                If subclustering increment cluster size by (default: 5).
  --maxsize MAXSIZE                    Maximum size of clusters for subclustering (default: 40).
  --heatmaps HEATMAPS                  Path to heatmap file (default: None).
  --pcutoff PCUTOFF                    Cutoff value for p-value if there is a p-value column in
                                       the heatmap.
  --squish LOWER,UPPER                 Squish the color scale to LOWER,UPPER bounds (default: None).
  --relwidths DEND,HEAT                Set relative widths of dendrogram and heatmap to DEND,HEAT
                                       (default: 1,1).
  --plotwidth PLOTWIDTH                Width of the dendrogram visualization (default: 30).



old shit
  --cut-threshold CUT_THRESHOLD, -t CUT_THRESHOLD
                        Cut the dendrogram at this threshold. Only used if
                        `--cut-method=hard`. (default: 0.3)
  --cut-method {dynamic,hard,none}, -m {dynamic,hard,none}
                        If `dynamic`, use dynamicTreeCut to determine clusters
                        without a hard threshold. If `none`, return all the
                        clusterings. Otherwise, use with `--cut-threshold` to
                        provide a numeric cut threshold. (default: dynamic)
  --dendrogram-style {rectangular,r,polar,p,none,n}, -s {rectangular,r,polar,p,none,n}
                        Plot the dendrogram in rectangular or polar
                        coordinates. If "none", then do not plot the dendrogram
                        (this is redundant with --no-plot). (default:
                        rectangular)
  --labels-use-clusters
  --labels-use-names [LABELS_USE_NAMES ...]
                        Label the dendrogram using columns from the nodetable.
                        This is a space-separated list of column names from the
                        nodetable. Pass columns as strings (column names).
                        (default: None)
  --labels-use-locs [LABELS_USE_LOCS ...]
                        Label the dendrogram using columns from the nodetable.
                        This is a space-separated list of integers indicating
                        columns from the nodetable (0-index, e.g., the first
                        column, which contains the node names, has index 0;
                        the second column has index 1, etc). (default: None)
  --labels-sep LABELS_SEP
                        The separator that will be used if multiple columns
                        from nodetable are used to label the dendrogram.
                        (default: | )
  --outdir OUTDIR       Save dendrogram and clusters to path. (default: None)
  --out-dendrogram OUT_DENDROGRAM, -d OUT_DENDROGRAM
                        Save dendrogram to path. (default: None)
  --no-plot             Do not plot the dendrogram. (default: False)
  --out-clusters OUT_CLUSTERS, -c OUT_CLUSTERS
                        Save clusters to path as tsv file with columns
                        "label", "cluster". When --threshold is 0 (the
                        default) each gene is put into a separate cluster
                        (i.e., every cluster has only a single gene).
                        (default: None)
  --no-clusters         Do not export clusters to file. (default: False)
  --path-to-conda-env PATH_TO_CONDA_ENV
  --path-to-rwrtoolkit PATH_TO_RWRTOOLKIT
  --multiplex MULTIPLEX
  --geneset GENESET
  --method METHOD
  --folds FOLDS
  --restart RESTART
  --tau TAU
  --numranked NUMRANKED
  --modname MODNAME
  --plot PLOT
  --threads THREADS
  --init-test-fullranks INIT_TEST_FULLRANKS
                        Create fullranks file for testing at the given path. (default: None)
  --verbose, -v         Default: WARNING; once: INFO; twice: DEBUG (default: 0)
  --version             Print version and exit. (default: False)
```

Examples
========

```sh
mentor \
    --distances /path/to/dissimilarity-matrix.tsv
    --outdir /path/to/outdir/
    --clusters CLUSTERS
    --map /path/to/map.txt
    --subcluster
    --increment INCREMENT
    --maxsize MAXSIZE
    --heatmaps /path/to/map.txt
    --squish=-1,1
    --relwidths=1,1
    --plotwidth 35
```

```sh
mentor \
    --rwr-fullranks /path/to/rwr-fullranks.tsv
    --partition
    --outdir /path/to/outdir/
    --clusters CLUSTERS
    --map /path/to/map.txt
    --subcluster
    --increment INCREMENT
    --maxsize MAXSIZE
    --heatmaps /path/to/map.txt
    --squish=-1,1
    --relwidths=1,1
    --plotwidth 35
```

```sh
mentor \
    --path-to-conda-env /path/to/conda/env
    --path-to-rwrtoolkit /path/to/RWRtoolkit
    --geneset /path/tp/geneset.txt
    --multiplex /path/to/multiplex.RData
    --partition
    --outdir /path/to/outdir/
    --clusters CLUSTERS
    --map /path/to/map.txt
    --subcluster
    --increment INCREMENT
    --maxsize MAXSIZE
    --heatmaps /path/to/map.txt
    --squish=-1,1
    --relwidths=1,1
    --plotwidth 35
```


If you have a gene set and a multiplex network, run RWR and partition the
geneset:

```sh
functional_partitioning \
    --path-to-rwrtoolkit <path to RWRtoolkit directory> \
    --multiplex <path to multiplex object> \
    --geneset <path to your gene set> \
    --outdir <path to output directory>
```

If you already have results from RWRtoolkit, you can pass in your 'fullranks'
file:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --outdir <path to output directory>
```

The default cut method is to use `dynamicTreeCut`, explicitly:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --cut-method dynamic \
    --outdir <path to output directory>
```

You can change this to use a hard threshold instead:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --cut-method hard \
    --cut-threshold 0.3 \
    --outdir <path to output directory>
```

You can change the labels on the leaves of the dendrogram. This requires that
you have a nodetable file. The nodetable is a TSV file with the first column as
the seed genes use for RWR-CV (singletons), e.g., something like this:

```
gene_id          gene_info
Potri.001G377800 This is gene 001G377800.
Potri.001G429430 This is gene 001G429430.
Potri.003G099600 This is gene 003G099600.
...
```

Then, provide the `--nodetable` and `--labels-use-{locs,names}` parameters to
tell this script which columns to use as labels on the dendrogram. E.g., to use
the column "gene_info" as the labels, you can either provide the column by name:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --nodetable <path to nodetable.tsv file> \
    --labels-use-names "gene_info"
```

Or provide the column by number:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --nodetable <path to nodetable.tsv file> \
    --labels-use-locs 1
```

If you want to keep the node ID in the label, provide that column as one of the
labels in addition to any others:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --nodetable <path to nodetable.tsv file> \
    --labels-use-names "gene_id" "gene_info"
```


Acknowledgements
================

Derived from the [cookiecutter data science][] project template.


<!-- LINKS -->

[cookiecutter data science]: https://drivendata.github.io/cookiecutter-data-science/
[software_release_life_cycle]: https://en.wikipedia.org/wiki/Software_release_life_cycle

<!-- END -->
