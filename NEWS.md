1.1.90
======

* Add options to generate violin plots instead of boxplots

* Add support for Jenks natural breaks and custom grouping to  step7TopGenesSurvival

* Add support for custom regions and MotifDb filtering options to step7LinkedDNAMethylationSitesMotifSearching

* step7TopGenesUserPeakOverlap now accepts peak names, generates results for each individual peak rather than each peak file, and accepts multiple input objects or directories

* Fix TENETSavedSizePlot bug that caused warnings and failure to respect an explicitly specified size

* Fix TENETCacheAllData to fail if any dataset does not exist on the hub

1.1.1
=====

* Merge survival bug fix from Bioconductor 3.21 release branch; see below

1.1.0
=====

* Release generated automatically by Bioconductor - no changes

1.0.1
=====

* Fix typo in survival calculation code which made the number of case samples with intermediate expression/methylation levels always 0

1.0.0
=====

* Initial Bioconductor release
