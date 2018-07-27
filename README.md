Replication Power
================
Nathan (Nat) Goodman
July 27, 2018

<!-- README.md is generated from README.Rmd. Please edit that file -->
*A collection of R scripts and documents exploring the power of replication to detect bad science. I treat replication as a statistical test, simulate proposed replication methods across a wide range of conditions, and estimate error rates for conditions of interest. The main point is that replication is a poor statistical test with unacceptable error rates under most conditions. It works as a validation test only when the original and replica studies are sampling nearly identical populations. Methods for testing whether the populations are similar work poorly under all conditions analyzed.*

**THE SOFTWARE IS STILL ROUGH and SOFTWARE DOCUMENTATION NONEXISTENT. PLEASE GET IN TOUCH IF YOU NEED HELP**

Overview
--------

The program explores the power of replication to detect bad science. The software simulates *studies* across a range of conditions, combines pairs of studies into *pairwise replications*, applies rules (called *measures*) for deciding which pairwise replications pass, summarizes the results as counts and pass rates, and finally computes false positive and negative rates for measures and conditions of interest. The main conclusion is that replication is a poor statistical test with unacceptable error rates under most conditions. Significance testing of the replica works fine as a validation test for exact and near-exact replications, but error rates increase rapidly as the populations diverge. All other tests have excessive error rates under all conditions analyzed. All tests have unacceptable error rates when used to check whether the original and replica studies are similar.

To calculate error rates, it's necessary to define explicit correctness criteria. The ones I use are

1.  *non-zero* - a replication instance is *true* if the population effect size of the first study is non-zero
2.  *same-effect* - a replication instance is *true* if both studies have the same population effect size; with *tolerance* *δ*, a replication instance is *true* if the two population effect sizes differ by at most *δ*

The measures appearing in this README document are

-   *sig2* - the replica study has a significant p-value
-   *sigm* - the fixed effect meta-analysis of the studies has a significant p-value
-   *d1.c2*, *d2.c1* - the standardized observed effect size (aka *Cohen's d*) of one study is in the confidence interval of the other
-   *d1.p2*, *d2.p1* - *Cohen's d* of one study is in the prediction interval of the other
-   *c1.c2* (resp. *p1.p2*) - the confidence (resp. prediction) intervals of the two studies overlap
-   *d2.scp1* - Uri Simonsohn's *small telescopes* method

All measures assume that the original study is significant (*sig1* in my notation) and the observed effect sizes of the two studies have the same sign. *Small telescopes* also assumes that *sig2* holds.

Installation and Usage
----------------------

The software is **not a package** and cannot be installed by `devtools::install_github` or related. Sorry. The simplest way to get the software is to download or clone the entire repo.

The code mostly uses base R capabilities but has a few dependencies: `RColorBrewer`, `akima`, and `pryr`. Since it's not a package, you have to manually install these packages if you don't already have them.

The recommended way to run the program is to `source` the file `R/repwr.R` into your R session; `R/repwr.R` will source the rest. Once loaded, you can run the program by executing the statement `run()` as shown below.

``` r
## This code block assumes your working directory is the root of the distribution.

source('R/repwr.R');
run();
```

This runs the program in a demo-like mode that quickly generates the simulated data and produces the figures that appear in this README document. The default computation simulates 625,000 replications and produces 16 figures. The simulation takes about 40 seconds on my small Linux server; the figures take about 60 seconds, much of which I think is spent rendering the plots over a remote X11 connection.

You can run each part separately by running one of the statements below.

``` r
## This code block assumes your working directory is the root of the distribution
## and you've already sourced R/repwr.R into your session

dodata();          # generate the simulated data
dodoc();           # generate the figures
dodoc(figsave=F);  # generate the figures without saving them
dodoc(fignew=F);   # generate and save the figures without plotting each in a new window
```

The program can also generate the data and figures for the other documents associated with the project: a blog post and (soon) a technical note with more details. **CAUTION: these take much longer to run**: about an hour each on my small Linux server. To generate these, execute `run()` with a suitable `doc` argument as shown below; you also need to specify `clean=T`, since by default the program reuses existing data for these large runs.

``` r
## This code block assumes your working directory is the root of the distribution.

source('R/repwr.R');
run(doc='repwr',clean=T);  # generate data and figures for blog post
run(doc='tech',clean=T);   # generate data and figures for technical note
```

Figures
-------

The default mode produces figures that illustrate the kinds of graphs the program can produce.

1.  line graphs showing error rates for a set of measures across chosen conditions - simple and intuitive but poor at showing data for too many measures and parameters
2.  heatmaps showing the same kind of data - still reasonably intuitive and somewhat better at depicting more measures and parameters
3.  rate-vs-rate scatter plots - able to display error rates across large swaths of parameter space but with less parameter resolution and perhaps less intuitive clarity
4.  aggregate line graphs - same data as rate-vs-rate scatter plots but for fewer measures and with better parameter resolution

The first group of figures are line graphs showing false positive and false negative rates for a few measures across a few conditions. The labels on the x-axis show the conditions: *n1*, *n2* are the sample sizes; *d1*, *d2* are the population effect sizes. The horizontal dashed lines demark the conventionally accepted thresholds of 0.05 for false positives and 0.20 for false negatives.

Figures 1-2 are for the *non-zero* correctness criterion; figures 3-4 are for *same-effect*.

<img src="figure/readme/m=1e3/figure_001_plotrate_nonzro_fpr.png" width="50%" /><img src="figure/readme/m=1e3/figure_002_plotrate_nonzro_fnr.png" width="50%" /><img src="figure/readme/m=1e3/figure_003_plotrate_sameff_fpr.png" width="50%" /><img src="figure/readme/m=1e3/figure_004_plotrate_sameff_fnr.png" width="50%" />

The next figures are heatmaps. Figures 5-6 show the same conditions as figures 1-2 but for more measures; figures 7-8 show more conditions. The red-to-blue transition is set at the conventionally accepted thresholds of 0.05 for false positives and 0.20 for false negatives. The dark vertical lines in figures 7-8 visually split each plot into separate "panels" for each value of *d2*.

<img src="figure/readme/m=1e3/figure_005_heatrate_nonzro_fpr.png" width="50%" /><img src="figure/readme/m=1e3/figure_006_heatrate_nonzro_fnr.png" width="50%" /><img src="figure/readme/m=1e3/figure_007_heatrate_nonzro_fpr_multi.png" width="50%" /><img src="figure/readme/m=1e3/figure_008_heatrate_nonzro_fnr_multi.png" width="50%" />

The next two figures (figures 9-10) are rate-vs-rate graphs for *exact* and *inexact* replications. Each point shows the mean false negative vs. mean false positive rate for specific conditions grouped by *n*1, *n*2. The dashed lines demark the conventionally acceptable error rates; the bottom left hand corner is the region where both error rates are acceptable. You'll note that for *exact*, *sig2* is the only measure with points in the acceptable region; for *inexact*, no points are in the acceptable region.

<img src="figure/readme/m=1e3/figure_009_roc_exact.png" width="50%" /><img src="figure/readme/m=1e3/figure_010_roc_inexact.png" width="50%" />

Figures 11-12 are aggregate line graphs showing the same data as the rate-vs-rate graphs above for fewer measures.

<img src="figure/readme/m=1e3/figure_011_rag_exact.png" width="50%" /><img src="figure/readme/m=1e3/figure_012_rag_inexact.png" width="50%" />

Recall that *sig2* works fine in exact replications but poorly in inexact ones (see figures 9-10). The next two figures (figures 13-14) show how *sig2* performs in *near exact* replications, ones where the population effect sizes differ slightly. The first is a rate-vs-rate graph showing *sig2* across various nearness values; the second is an aggregate line graph showing the same data.

<img src="figure/readme/m=1e3/figure_013_multi_sig2_rocm.png" width="50%" /><img src="figure/readme/m=1e3/figure_014_multi_sig2_ragm.png" width="50%" />

The final two figures (figures 15-16) compare *sig2* and *d2.scp1* (Uri Simonsohn's *small telescopes* method). The differences are quite small.

<img src="figure/readme/m=1e3/figure_015_small_telescopes_roc.png" width="50%" /><img src="figure/readme/m=1e3/figure_016_small_telescopes_rag.png" width="50%" />

See Also
--------

A blog post discussing the approach and results is available in [html](https://natgoodman.github.io/repwr/repwr.stable.html) and [pdf](https://natgoodman.github.io/repwr/repwr.stable.pdf) on the [GitHub Pages site](https://natgoodman.github.io/repwr) associated with this repository and will soon be posted on a blog site TBD. It's also in the repository as files [repwr.stable.html](https://github.com/natgoodman/repwr/repwr.stable.html) and [repwr.stable.pdf](https://github.com/natgoodman/repwr/repwr.stable.pdf). (But note that GitHub, unlike GitHub Pages, renders html files as raw text).

A document with technical details will soon be available in [html](https://natgoodman.github.io/repwr/tech.stable.html) and [pdf](https://natgoodman.github.io/repwr/tech.stable.pdf) on the [GitHub Pages site](https://natgoodman.github.io/repwr) and in the repository as files [tech.stable.html](https://github.com/natgoodman/repwr/tech.stable.html) and [tech.stable.pdf](https://github.com/natgoodman/repwr/tech.stable.pdf).

Author
------

Nathan (Nat) Goodman, (natg at shore.net)

Bugs and Caveats
----------------

Please report any bugs, other problems, and feature requests using the [GitHub Issue Tracker](https://github.com/natgoodman/repwr/issues). I will be notified, and you'll be apprised of progress. As already noted, the software is still rough and software documentation nonexistent.

Copyright & License
-------------------

Copyright (c) 2018 Nathan Goodman

The software is **open source and free**, released under the [MIT License](https://opensource.org/licenses/MIT). The documentation is **open access**, released under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0).
