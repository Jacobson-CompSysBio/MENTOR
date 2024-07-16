**MENTOR**: **M**echanistic **E**xploration of **N**etworks for **T**eam-based **O**mics **R**esearch
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
usage: mentor [-h] [--geneset /path/to/geneset.txt] [--multiplex /path/to/multiplex.RData] [--outdir /path/to/outdir]
                   [--threads threads] [--verbose] [--version] [--distances /path/to/dissimilarity-matrix.tsv] [--clusters clusters]
                   [--map /path/to/map.txt] [--subcluster] [--increment increment] [--maxsize maxsize] [--heatmaps /path/to/heatmap.txt]
                   [--reordercols] [--pcutoff pcutoff] [--squish=lower,upper] [--relwidths=dend,heat] [--plotwidth plotwidth]

arguments:

  -h, --help                                          Show this help message and exit
  --outdir /path/to/outdir                            Save output to path (default: None)
  --outfile outfile                                   Description to append to filenames (default: None)
  --multiplex /path/to/multiplex.RData                Path to multiplex network (default: None)
  --geneset /path/to/geneset.txt                      Path to gene set file (default: None)
  --threads threads                                   Number of threads to use in calculation (default: based on system) 
  --verbose, -v                                       Gives output from runs. Default: WARNING; once: INFO; twice: DEBUG (default: 0)
  --version                                           Print version and exit (default: False)
  --distances /path/to/dissimilarity-matrix.tsv       Path to dissimilarity-matrix.tsv (default: None)
  --clusters clusters                                 Number of clusters for dendrogram (default: 3)
  --map /path/to/map.txt                              Path to gene symbol to gene ID mapping file (default: None)
  --subcluster                                        Subcluster the dendrogram (default: False)   
  --increment increment                               Increment subclustering cluster size by this value (default: 5)
  --maxsize maxsize                                   Maximum size of clusters for subclustering (default: 40)
  --heatmaps /path/to/heatmap.txt                     Path to heatmap file if including heatmaps for log2FC or effect size (default: None)
  --plottype shape                                    Type of dendrogram (either rectangular or polar; default: rectangular)
  --reordercols                                       Reorder columns of heatmap with clustering (default: False)
  --legendtitle legendtitle                           Titles to give to legends; for rectangular dendrogram title to give to continuous legend;
                                                      for polar dendrogram titles to give to continuous legend and factor legend (default: log2FC)
  --squish lower,upper                                Squish the color scale to LOWER,UPPER bounds (default: None)
  --pcutoff pcutoff                                   Cutoff value for p-value if there is a p-value column in
                                                      the heatmap
  --relwidths dend,heat                               Set relative widths of dendrogram and heatmap to dend,heat
                                                      (default: 1,1)
  --clusterlabelsize labelsize                        Size of cluster labels associated with polar dendrogram (default: 2.5)
  --heatmaplabelsize heatmaplabelsize                 Size of the labels associated with the polar dendrogram heatmap (gene names; default: 0.75)
  --groupcolors #00000F,#FFFFFF,#000FFF               Colors for the group blocks in polar dendrogram (default: None)
  --trackheight  height,height,height                 Width of the tracks in the polar dendrogram (heatmap1,heatmap2,dendrogram; default: 0.2, 0.2, 0.2)
  --highlightindex 0,1,2                              Sector(s) to highlight on polar dendrogram (sector1,sector2,...,sectorN; default: None)
  --highlightcolor #34EBDC                            Color to use for highlighted sectors on polar dendrogram (default: #34EBDC)
  --plotwidth plotwidth                               Width of the dendrogram visualization (default: 30)
  --plotheight plotheight                             Height of the dendrogram/heatmap visualization (default: None)
```

Examples
========

The `geneset.txt` and `multiplex.RData` are required to run `mentor` to obtain a dissimilarity matrix and dendrogram visualization. The `map.txt` and `heatmap.txt` are not required to run `mentor`.

The `geneset.txt` table should be a tab-separated text file with **no header**. Three columns should be present for the project name *(string)*, ensembl IDs *(string)*, and weights *(numeric)*. The example below displays a geneset table with five genes. 

|              |                 |   |
| ------------ | --------------- | - |
| project_name | ENSG00000008311 | 1 |
| project_name | ENSG00000141338 | 1 |
| project_name | ENSG00000172350 | 1 |
| project_name | ENSG00000042980 | 1 |
| project_name | ENSG00000198099 | 1 |

The `map.txt` table can be used to assign different labels to the dendrogram branches. The table should also be a tab-separated text file with a **header**. Two columns should be present for the ensembl ID *(string)* and associated label *(string)* that you would like to display in the dendrogram branch labels. The example below displays a mapping table for the same five genes. If you would like to include a heatmap in the visualization then you must ensure that the labels in the map table match the labels in the heatmap table.

|     ensembl     |  label  |
| --------------- | ------- |
| ENSG00000008311 |   AASS  |
| ENSG00000141338 |  ABCA8  |
| ENSG00000172350 |  ABCG4  |
| ENSG00000042980 | ADAM28  |
| ENSG00000198099 |   ADH4  |

The `heatmap.txt` table is required if the user would like to include a heatmap next to the dendrogram. The table should also be a tab-separated text file with a **header**. Three columns should be present for the label *(string)*, value *(numeric)*, and data source *(string)*. Each unique data source will be presented as a new column in the heatmap. The order of the source columns from left to right in the heatmap will be the same order as their order in the `heatmap.txt` table from top to bottom. If you would like to reorder the columns by similarity then specify `--reordercols`.The total number of columns is dependent on the type of information that you would like to present in the heatmap. The example below displays a heatmap table for the same five genes where we have bulk RNA-seq and GWAS data sources associated with these genes. You can see that all five genes were implicated in the RNA-seq data source but only three were implicated in the GWAS data source.

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

To run `mentor` on a gene set and multiplex network:

```sh
mentor \
    --geneset </path/tp/geneset.txt> \
    --multiplex </path/to/multiplex.RData> \
    --outdir </path/to/outdir/>
```

To run `mentor` on a gene set and multiplex and have custom dendrogram labels:

```sh
mentor \
    --geneset </path/tp/geneset.txt> \
    --multiplex </path/to/multiplex.RData> \
    --outdir </path/to/outdir/> \
    --map </path/to/map.txt/>
```

To run `mentor` on a gene set and multiplex and have custom dendrogram labels and a heatmap:

```sh
mentor \
    --geneset </path/tp/geneset.txt> \
    --multiplex </path/to/multiplex.RData> \
    --outdir </path/to/outdir/> \
    --map </path/to/map.txt/> \
    --heatmaps </path/to/heatmap.txt/>
```

To customize a dendrogram using a `dissimilarity-matrix.tsv`: 

```sh
mentor \
    --distances </path/to/dissimilarity-matrix.tsv> \
    --outdir </path/to/outdir/> \
    --clusters <clusters> \
    --map </path/to/map.txt> \
    --subcluster \
    --increment <increment> \
    --maxsize <maxsize> \
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
