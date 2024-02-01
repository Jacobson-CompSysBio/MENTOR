MENTOR: Mechanistic Exploration of Networks for Team-based Omics Research
=======================

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

usage: mentor [-h] [--geneset GENESET] [--multiplex MULTIPLEX] [--outdir OUTDIR]
                   [--threads THREADS] [--verbose] [--version] [--distances DISTANCES] [--clusters CLUSTERS]
                   [--map MAP] [--subcluster] [--increment INCREMENT] [--maxsize MAXSIZE] [--heatmaps HEATMAPS]
                   [--pcutoff PCUTOFF] [--squish LOWER,UPPER] [--relwidths DEND,HEAT] [--plotwidth PLOTWIDTH]

arguments:
  -h, --help                           Show this help message and exit.
  --geneset GENESET                    Path to gene set file.
  --multiplex MULTIPLEX                Path to multiplex network.
  --outdir OUTDIR                      Save output to path (default: None).
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

The `geneset.txt` table should be a tab-separated text file with **no header**. Three columns should be present for the project name *(string)*, ensembl IDs *(string)*, and weights *(numeric)*. The example below displays a geneset table with five genes. 

|              |                 |   |
| ------------ | --------------- | - |
| project_name | ENSG00000008311 | 1 |
| project_name | ENSG00000141338 | 1 |
| project_name | ENSG00000172350 | 1 |
| project_name | ENSG00000042980 | 1 |
| project_name | ENSG00000198099 | 1 |

The `map.txt` table should also be a tab-separated text file with a **header**. Two columns should be present for the ensembl ID *(string)* and associated label *(string)* that you would like to display in the dendrogram branch labels. The example below displays a mapping table for the same five genes. If you would like to include a heatmap in the visualization then you must ensure that the labels in the map table match the labels in the heatmap table.

|     ensembl     |  label  |
| --------------- | ------- |
| ENSG00000008311 |   AASS  |
| ENSG00000141338 |  ABCA8  |
| ENSG00000172350 |  ABCG4  |
| ENSG00000042980 | ADAM28  |
| ENSG00000198099 |   ADH4  |

The `heatmap.txt` table should also be a tab-separated text file with a **header**. Three columns should be present for the label *(string)*, value *(numeric)*, and data source *(string)*. Each unique data source will be presented as a new column in the heatmap. The total number of columns is dependent on the type of information that you would like to present in the heatmap. The example below displays a heatmap table for the same five genes where we have bulk RNA-seq and GWAS data sources associated with these genes. You can see that all five genes were implicated in the RNA-seq data source but only three were implicated in the GWAS data source.

|  label  |  value  |  source  |
| ------- | ------- | -------- |
|   AASS  |   1.5   |  RNA-seq |
|  ABCA8  |   2.5   |  RNA-seq |
|  ABCG4  |   0.3   |  RNA-seq |
| ADAM28  |   -1.1  |  RNA-seq |
|   ADH4  |   -0.8  |  RNA-seq |
|   AASS  |    1    |   GWAS   |
|  ABCA8  |    1    |   GWAS   |
|   ADH4  |    1    |   GWAS   |

To run RWRtoolkit on a gene set and multiplex network followed by mentor:

```sh
mentor \
    --geneset </path/tp/geneset.txt> \
    --multiplex </path/to/multiplex.RData> \
    --outdir </path/to/outdir/> \
    --clusters <CLUSTERS> \
    --map </path/to/map.txt> \
    --subcluster \
    --increment <INCREMEN> \
    --maxsize <MAXSIZE> \
    --heatmaps </path/to/heatmap.txt> \
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
    --heatmaps </path/to/heatmap.txt> \
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
