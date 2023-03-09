% Functional Partitioning
% J. Izaak Miller

Functional Partitioning
=======================

- Integrate multi-omics data using multiplex networks.
- Identify functionally-related groups of genes using random walk with restart
  on multiplex networks.

Installation
============

```sh
$ git clone https://github.com/izaakm/jail-functional-partitioning
$ cd jail-functional-partitioning
$ pip install .
```

Usage
=====

```
usage: functional_partitioning [-h] [--rwr-fullranks RWR_FULLRANKS]
                               [--nodetable NODETABLE] [--partition]
                               [--no-partition]
                               [--cut-threshold CUT_THRESHOLD]
                               [--cut-method {dynamic,hard,none}]
                               [--dendrogram-style {rectangular,r,polar,p,none,n}]
                               [--no-plot] [--labels-use-clusters]
                               [--labels-use-names [LABELS_USE_NAMES ...]]
                               [--labels-use-locs [LABELS_USE_LOCS ...]]
                               [--labels-sep LABELS_SEP] [--outdir OUTDIR]
                               [--out-dendrogram OUT_DENDROGRAM]
                               [--out-clusters OUT_CLUSTERS]
                               [--path-to-conda-env PATH_TO_CONDA_ENV]
                               [--path-to-rwrtoolkit PATH_TO_RWRTOOLKIT]
                               [--multiplex MULTIPLEX] [--geneset GENESET]
                               [--method METHOD] [--folds FOLDS]
                               [--restart RESTART] [--tau TAU]
                               [--numranked NUMRANKED] [--modname MODNAME]
                               [--plot PLOT] [--threads THREADS]
                               [--init-test-fullranks INIT_TEST_FULLRANKS]
                               [--verbose] [--version]

Partition seeds from `RWR-CV --method=singletons ...` into clusters.

options:
  -h, --help            show this help message and exit
  --rwr-fullranks RWR_FULLRANKS, -f RWR_FULLRANKS
                        Path to "fullranks" file from `RWR-CV
                        --method=singletons ...` (default: None)
  --nodetable NODETABLE
                        Path to "nodetable" file. This is a TSV file where the
                        first column is the node name (i.e., the seed genes
                        from RWR-fullranks). (default: None)
  --partition, -p       Perform functional partitioning on "seed genes" from
                        RWR fullranks file. This is the default. (default:
                        True)
  --no-partition        Do not perform functional partitioning. (default:
                        True)
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
                        coordinates. If "none", then do not plot the
                        dendrogram (this is redundant with --no-plot).
                        (default: rectangular)
  --no-plot             Do not plot the dendrogram. (default: False)
  --labels-use-clusters
  --labels-use-names [LABELS_USE_NAMES ...]
                        Label the dendrogram using columns from the nodetable.
                        This is a space-separated list of column names from
                        the nodetable. Pass columns as strings (column names).
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
  --out-clusters OUT_CLUSTERS, -c OUT_CLUSTERS
                        Save clusters to path as tsv file with columns
                        "label", "cluster". When --threshold is 0 (the
                        default) each gene is put into a separate cluster
                        (i.e., every cluster has only a single gene).
                        (default: None)
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
                        Create fullranks file for testing at the given path.
                        (default: None)
  --verbose, -v         Default: WARNING; once: INFO; twice: DEBUG (default:
                        0)
  --version             Print version and exit. (default: False)
```

Examples
========

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
    --partition \
    --nodetable <path to nodetable.tsv file> \
    --labels-use-names "gene_info"
```

Or provide the column by number:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --partition \
    --nodetable <path to nodetable.tsv file> \
    --labels-use-locs 1
```

If you want to keep the node ID in the label, provide that column as one of the
labels in addition to any others:

```sh
functional_partitioning \
    --rwr-fullranks <path to fullranks file> \
    --partition \
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
