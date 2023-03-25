% Changelog

<!-- # [<version>] - <date> -->
<!-- [<version>]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v<version> -->
<!-- ## Added for new features. -->
<!-- ## Changed for changes in existing functionality. -->
<!-- ## Deprecated for soon-to-be removed features. -->
<!-- ## Fixed for any bug fixes. -->
<!-- ## Removed for now removed features. -->
<!-- ## Security in case of vulnerabilities. -->

# [Unreleased]

<!-- # [<version>] - <date> -->
<!-- [<version>]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v<version> -->
<!-- ## Added -->
<!-- ## Changed -->
<!-- ## Deprecated -->
<!-- ## Fixed -->
<!-- ## Removed -->
<!-- ## Security -->


# [0.6.0] - 2023/03/25


## Added

- Add basic summary statistics for pairwise distances
    - Calculate basic summary statistics for pairwise distances, dump to file
- Add `__all__` to `__init__.py`
- Add `_version.py`, `_plot.py`
- Add tests for CLI entry point, etc:
    `test_basic.py`,
    `test_cli.py`,
    `test_clustering.py`,
    `test_rwrtoolkit.py`,
    `test_sample_data.py`
- Split up `functional_partitioning` module:
    `_cluster.py`,
    `_datasets.py`,
    `_metrics.py`,
    `_plot.py`,
    `_rwrtoolkit.py`,
    `_utils.py`,
    `_version.py`,
    `functional_partitioning.py` -> `cli.py`
- Makefile
- `get_scores_vs_ranks_curve` function for use with `get_elbow`
- `--no-clusters` option
- Examples:
    - Add example of `--cut-method` and `--cut-threshold`
    - Add example using RWRtoolkit


## Changed

- Allow passing `--cut-method=none` to return all clusterings
- Clarify setting `out_{clusters,dendrogram}` variables
- Move `calc_chi` function to `_metrics.py`
- Move `calc_threshold` function to `_cluster.py`
- Move `get_elbow` function to `_metrics.py`
- Move plotting functions to `_plot.py`
- Rename `functional_partitioning.py` to `cli.py`
- Rename `get_top_ranked` function
- Restructure package: prefix modules with underscore
- Replace `--threshold` with `--cut-threshold` and `--cut-method`
- Use `__version__` as string instead of `get_version`
- Wrap `_plot_dendrogram_{polar,rectangular}` w/ `draw_dendrogram`
- Create the output directory (but not its parents)
- Attempt to clarify installation instructions
- Refactor tests (use class `unittest.TestCase`)
- Include `dynamicTreeCut` in environment.yml for conda installations
- Rank transformation of scores: use `method='first'`
- Require dissimilarities as a condensed distance matrix
- Adjust dendrogram figsize depending on orientation
- Pass optimal_ordering to scipy's `hierarchy.linkage` fxn
- Set `min_cluster_size_` to 1 instead of 3
    - Dynamic tree cut tries to assign genes that really should be singletons to clusters. This minimizes that behavior, but it's still a problem.
- Use `__version__` as string instead of `get_version`
- Wrap `_plot_dendrogram_{polar,rectangular}` w/ `draw_dendrogram`
- Clarify setting `out_{clusters,dendrogram}` variables
- Move `calc_threshold` fxn to `cluster`
- Move `calc_chi` fxn to `metrics`
- Move `get_elbow` fxn to `metrics`
- Use the `HierarchicalClustering` class with `dynamicTreeCut` (cutreeHybrid) as the default clustering method.


## Removed

- Remove `get_top_ranked` fxn


## Fixed

- Get version from `_version.__version__`
- Install `seaborn`
- `import dynamicTreeCut from dynamicTreeCut`
- Adhere to standard test name conventions
- Restructure package: prefix modules with underscore to avoid 'flat-layout' error:
  ```
  error: Multiple top-level modules discovered in a flat-layout:
  ['metrics', 'datasets', 'cluster', 'functional_partitioning', 'rwrtoolkit'].
  ```
    - Import modules with leading underscores
- Calculate the elbow using mean of scores-vs-ranks vectors
    - Bug in elbow calculation introduced during recent refactoring.
    - Sort mean scores before calculating elbow
- Modify package structure for install
- Write clusters to clusters.tsv file


# [0.5.1] - 2023/03/08

## Added

- `metric` parameter was added in sklearn `AgglomerativeClustering` (upstream).
- `dynamicTreeCut` dependency for `HierarchicalClustering`
- `check_symmetry` function to verify square (uncondensed) distance matrix

## Changed

- Require `scikit-learn>=1.2.0`
- `HierarchicalClustering` : Use `scipy.cluster.hierarchy.{linkage,cut_tree}` to create `labels_`.
- `HierarchicalClustering` : Use `metric` parameter instead of `affinity`.
- `HierarchicalClustering` : Use `dynamicTreeCut.cutreeHybrid` as default method for creating clusters.

## Deprecated

- `affinity` parameter was decrepated in sklearn `AgglomerativeClustering` (upstream).

<!-- ## Removed for now removed features. -->
<!-- ## Fixed for any bug fixes. -->
<!-- ## Security in case of vulnerabilities. -->

# [0.4.0] - 2023/01/09

Initial release.

## Added

- CHANGELOG.md
- `cluster.HierarchicalClustering`
- RWRtoolkit wrapper functions (`rwrtoolkit.rwr_singletons`)

## Changed

- Split loading 'fullranks' and converting it to 'X' matrix into separate
  functions

---

See also:

- [About releases] 
- [Example changelog]


<!-- LINKS -->

[unreleased]: https://github.com/izaakm/jail-functional-partitioning/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v0.6.0
[0.5.1]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v0.5.1
[0.4.0]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v0.4.0
[About releases]: https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases
[Example changelog]: https://github.com/Level/level/blob/master/CHANGELOG.md

<!-- END -->
