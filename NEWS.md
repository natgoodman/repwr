Revision history for repwr repository
================
Nathan (Nat) Goodman
May 3, 2018

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->
Release 0.93 2018-05-dd
-----------------------

-   Removed obsolete scripts including scope.R and distr.R
-   Added current script: repwr.R
-   Split repwr.R into separate files
-   Removed other obsolete files
-   For data, added current top-level files needed for analysis and plotting. Did not add data subdirectories as these are huge (261MB for m=1e3, 980MB for m=1e4) and not needed once the simulation is done

Release 0.92 2018-01-01
-----------------------

### scope.R

-   Added plot\_cross, plot\_pclose1
-   Integrated d2t probability functions from distr.R
-   Systematically precomputed useful data tables: d1, d33, dcross, dsig, dclose, dcrit

### distr.R

-   Added pval2t, pval2d functions
-   Corrected usage of lower.tail parameter in d2t probability functions

Release 0.91 2017-12-19
-----------------------

-   Change sampling distribution from normal to noncentral t.
-   Remove unused code.
-   Clean up code for translating Cohen's d to t-statistic and t-distribution; simplify function names.

Release 0.90 2017-12-07
-----------------------

First version. Ready for feedback from external readers.

Copyright & License
-------------------

Copyright (c) 2017-2018 Nathan Goodman

The software is **open source and free**, released under the [MIT License](https://opensource.org/licenses/MIT). The documentation is **open access**, released under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0).
