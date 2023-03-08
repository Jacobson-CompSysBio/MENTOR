% Changelog

<!-- # [<version>] - <date> -->
<!-- ## Added for new features. -->
<!-- ## Changed for changes in existing functionality. -->
<!-- ## Deprecated for soon-to-be removed features. -->
<!-- ## Removed for now removed features. -->
<!-- ## Fixed for any bug fixes. -->
<!-- ## Security in case of vulnerabilities. -->
<!-- [<version>]: https://github.com/Level/level/releases/tag/<version> -->

# [Unreleased]

...

# [0.5.0] - 2023/03/08

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

[unreleased]: https://github.com/izaakm/jail-functional-partitioning/compare/v0.5.1...HEAD
[0.5.1]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v0.5.1
[0.4.0]: https://github.com/izaakm/jail-functional-partitioning/releases/tag/v0.4.0
[About releases]: https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases
[Example changelog]: https://github.com/Level/level/blob/master/CHANGELOG.md

<!-- END -->
