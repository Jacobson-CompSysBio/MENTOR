MENTOR: Mechanistic Exploration of Networks for Team-based Omics Research
=======================

- Integrate multi-omics data using multiplex networks
- Identify functionally-related groups of genes using random walk with restart
  on multiplex networks

Installation
============

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

The geneset table should be a tab-separated text file with *no header*. Three columns should be present for the project name (string), ensembl IDs (string), and weights (numeric). The example below displays a geneset table with five genes. 

|              |                 |   |
| ------------ | --------------- | - |
| project_name | ENSG00000008311 | 1 |
| project_name | ENSG00000141338 | 1 |
| project_name | ENSG00000172350 | 1 |
| project_name | ENSG00000042980 | 1 |
| project_name | ENSG00000198099 | 1 |

The map table should also be a tab-separated text file with a header. Two columns should be present for the ensembl ID (string) and associated label (string) that you would like to display in the dendrogram branch labels. The example below displays a mapping table for the same five genes. If you would like to include a heatmap in the visualization then you must ensure that the labels in the map table match the labels in the heatmap table.

|     ensembl     |  label  |
| --------------- | ------- |
| ENSG00000008311 |   AASS  |
| ENSG00000141338 |  ABCA8  |
| ENSG00000172350 |  ABCG4  |
| ENSG00000042980 | ADAM28  |
| ENSG00000198099 |   ADH4  |

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



Acknowledgements
================

Derived from the [cookiecutter data science][] project template.


<!-- LINKS -->

[cookiecutter data science]: https://drivendata.github.io/cookiecutter-data-science/
[software_release_life_cycle]: https://en.wikipedia.org/wiki/Software_release_life_cycle

<!-- END -->
