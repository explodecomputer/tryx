*Note: Apologies if the scan function times out sometimes - there are changes being made to the underlying database which will be resolved soon. We will also be releasing a major update to this package in the next few days that will make it more resilient to timeouts*


# Treasure your exceptions! (TRYX)

<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![Travis-CI build status](https://travis-ci.org/explodecomputer/tryx.svg?branch=master)](https://travis-ci.org/explodecomputer/tryx) [![codecov](https://codecov.io/github/explodecomputer/tryx/branch/master/graphs/badge.svg)](https://codecov.io/github/explodecomputer/tryx)

<!-- badges: end -->

**Major update:**

This package has been under major development and the way in which it is implemented has changed quite substantially. See below for how to revert to the previous version if necessary


---

**Full documentation including a vignette is available here:** https://explodecomputer.github.io/tryx/

---

This package will perform MR-TRYX analysis, which entails the following. 

In MR analysis a major assumption is that the SNP influences the outcome only through the exposure. If the SNP influences the outcome through some other (candidate) traits, in addition to influencing through the exposure, then the MR estimate can be biased.

Knowing which traits through which SNPs might be acting in a horizontal pleiotropic manner allows us to exploit these outliers to:

1. Identify novel candidate traits that influence the outcome (and possibly the exposure also)
2. Adjust the SNP-outcome / SNP-exposure ratio based on knowledge of the alternative pathways, thereby reducing heterogeneity in the original exposure-outcome effect estimate.

This package uses the MR-Base infrastructure to search for alternative pathways through which instruments might influence the exposure.

---

The package name TRYX (pronounced 'tricks') is taken from the phrase 'TReasure Your Exceptions', a quote from William Bateson (1908). 

## Requirements

This software is an R package that depends on various other packages available on CRAN (see DESCRIPTION file). It has only been tested on some R versions >= 3.2.0. The software will run on any standard laptop, desktop or server for which R can be installed.

## Installation

A beta development version can be installed directly from this repository.

First install the `TwoSampleMR` and `RadialMR` R packages:

```r
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("WSpiller/RadialMR")
```

Next install the `tryx` package:

```r
devtools::install_github("explodecomputer/tryx")
```

You may also want to install some plotting packages

```r
install.packages(c("ggplot2", "ggrepel", "igraph"))
```

and a package for simulating genotype-phenotype maps

```r
devtools::install_github("explodecomputer/simulateGP")
```

It should not take more than a few minutes to install all of these packages.

### Installing previous versions

You can go back to an earlier version using:

```r
devtools::install_github("explodecomputer/tryx@0.1.1")
```

## Citation

If you using MR-TRYX R package:

[Cho Y, Haycock P, Sanderson E, Gaunt T, Zheng J, Morris A, Davey Smith G, Hemani G. </br>
**MR-TRYX: A Mendelian randomization framework that exploits horizontal pleiotropy to infer novel causal pathways.** <br/>
Nature communications [Accepted]. Current version is available at on BioRxiv.](https://www.biorxiv.org/content/10.1101/476085v3)

For the IEU GWAS database, MR-Base or the TwoSamleMR R package:
[Hemani G, Zheng J, Elsworth B, Wade KH, Baird D, Haberland V, Laurin C, Burgess S, Bowden J, Langdon R, Tan VY, Yarmolinsky J, Shihab HA, Timpson NJ, Evans DM, Relton C, Martin RM, Davey Smith G, Gaunt TR, Haycock PC, The MR-Base Collaboration.</br>
**The MR-Base platform supports systematic causal inference across the human phenome.** <br/>
eLife 2018;7:e34408. doi: 10.7554/eLife.34408](https://elifesciences.org/articles/34408)



