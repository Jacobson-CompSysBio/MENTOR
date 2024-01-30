MENTOR: Mechanistic Exploration of Networks for Team-based Omics Research
=======================

- Integrate multi-omics data using multiplex networks.
- Identify functionally-related groups of genes using random walk with restart
  on multiplex networks.

Installation
============

Standard installation:

```sh
$ git clone <link>
$ cd MENTOR
$ make
```

If you prefer to install dependencies with conda, you can use the provided
`environment.yml` file:

```sh
$ git clone <link>
$ cd MENTOR
$ conda env update -f ./environment.yml
$ conda activate mentor
$ make
```

Or create a new environment called 'mentor':

```sh
$ git clone <link>
$ cd MENTOR
$ conda env create -f ./environment.yml
$ conda activate mentor
$ make
```

Verify the installation:

```sh
$ mentor --help
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
```

Examples
========

To run RWRtoolkit on a gene set and multiplex network followed by mentor:

```sh
mentor \
    --geneset </path/tp/geneset.txt> \
    --multiplex </path/to/multiplex.RData> \
    --partition \
    --outdir </path/to/outdir/> \
    --clusters <CLUSTERS> \
    --map </path/to/map.txt> \
    --subcluster \
    --increment <INCREMEN> \
    --maxsize <MAXSIZE> \
    --heatmaps </path/to/map.txt> \
    --squish=<-1,1> \
    --relwidths=<1,1> \
    --plotwidth <35> \
```

To customize a dendrogram from a mentor dissimilarity-matrix: 

```sh
mentor \
    --distances </path/to/dissimilarity-matrix.tsv> \
    --outdir </path/to/outdir/> \
    --clusters <CLUSTERS> \
    --map </path/to/map.txt> \
    --subcluster \
    --increment <INCREMEN> \
    --maxsize <MAXSIZE> \
    --heatmaps </path/to/map.txt> \
    --squish=<-1,1> \
    --relwidths=<1,1> \
    --plotwidth <35> \
```




old shit

If you have a gene set and a multiplex network, run RWR and mentor the geneset:

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
