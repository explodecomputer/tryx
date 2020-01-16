# Treasure your exceptions! (TRYX)

This package will perform MR-TRYX analysis, which entails the following. 

In MR analysis a major assumption is that the SNP influences the outcome only through the exposure. If the SNP influences the outcome through some other (candidate) traits, in addition to influencing through the exposure, then the MR estimate can be biased.

Knowing which traits through which SNPs might be acting in a horizontal pleiotropic manner allows us to exploit these outliers to:

1. Identify novel candidate traits that influence the outcome (and possibly the exposure also)
2. Adjust the SNP-outcome / SNP-exposure ratio based on knowledge of the alternative pathways, thereby reducing heterogeneity in the original exposure-outcome effect estimate.

This package uses the MR-Base infrastructure to search for alternative pathways through which instruments might influence the exposure.

---

The package name TRYX (pronounced 'tricks') is taken from the phrase 'TReasure Your Exceptions', a quote from William Bateson (1908). 

--- 

## Requirements

This software is an R package that depends on various other packages available on CRAN (see DESCRIPTION file). It has only been tested on some R versions >= 3.2.0. The software will run on any standard laptop, desktop or server for which R can be installed.

## Installation

A beta development version can bbe installed directly from this repository.

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

## Citation

