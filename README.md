Functional Partitioning
==============================

Functional partitioning.

Installation
------------

```sh
$ git clone https://github.com/izaakm/jail-functional-partitioning
$ cd jail-functional-partitioning
$ pip install .
```

Usage
-----

```
usage: functional_partitioning.py
    [-h] [--rwr-fullranks RWR_FULLRANKS] [--partition] [--threshold THRESHOLD]
    [--dendrogram-style {rectangular,r,polar,p}] [--out-dir OUT_DIR]
    [--out-dendrogram OUT_DENDROGRAM] [--out-clusters OUT_CLUSTERS] [--verbose]

Partition seeds from `RWR-CV --method=singletons ...` into clusters.

options:
    -h, --help
        show this help message and exit
    --rwr-fullranks RWR_FULLRANKS, -f RWR_FULLRANKS
        Path to "fullranks" file from `RWR-CV --method=singletons ...`
        (default: None)
    --partition, -p
        [PLACEHOLDER] Perform functional partitioning on "seed genes" from RWR
        fullranks file. This is the default. (default: True)
    --threshold THRESHOLD, -t THRESHOLD
        Apply threshold to dendrogram. Genes in branches below this threshold
        will be grouped into clusters; other genes are considered isolates
        (separate clusters, each with a single gene). Value can be float or
        "mean". If the value is "mean", then use the mean branch height as the
        cluster threshold; this can be useful for a first pass. (default: 0)
    --dendrogram-style {rectangular,r,polar,p}, -s {rectangular,r,polar,p}
        Plot the dendrogram in rectangular or polar coordinates. Default is
        rectangular. (default: rectangular)
    --out-dir OUT_DIR
        Save dendrogram and clusters to path. (default: None)
    --out-dendrogram OUT_DENDROGRAM, -d OUT_DENDROGRAM
        Save dendrogram to path. (default: None)
    --out-clusters OUT_CLUSTERS, -c OUT_CLUSTERS
        Save clusters to path as tsv file with columns "label", "cluster". When
        --threshold is 0 (the default) each gene is put into a separate cluster
        (i.e., every cluster has only a single gene). (default: None)
    --verbose, -v
        Default: WARNING; once: INFO; twice: DEBUG (default: 0)
```



Acknowledgements
----------------

Derived from the [cookiecutter data science][] project template.


<!-- LINKS -->

[cookiecutter data science]: https://drivendata.github.io/cookiecutter-data-science/

<!-- END -->
