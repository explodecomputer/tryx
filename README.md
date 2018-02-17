# Treasure your exceptions (TRYX)

This package will perform MR-TRYX analysis, which entails the following. 

In MR analysis a major assumption is that the SNP influences the outcome only through the exposure. If the SNP influences the outcome through some other (candidate) traits, in addition to influencing through the exposure, then the MR estimate can be biased.

Knowing which traits through which SNPs might be acting in a horizontal pleiotropic manner allows us to exploit these outliers to:

1. Identify novel candidate traits that influence the outcome (and possibly the exposure also)
2. Adjust the SNP-outcome / SNP-exposure ratio based on knowledge of the alternative pathways, thereby reducing heterogeneity in the original exposure-outcome effect estimate.

This package uses the MR-Base infrastructure to search for alternative pathways through which instruments might influence the exposure.

---

The package name TRYX (pronounced 'tricks') is taken from the phrase 'TReasure Your Exceptions', a quote from William Bateson. 

--- 

## Installation

First install the `TwoSampleMR` and `RadialMR` R package:

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

---

## Guide

### Basic analysis

Begin by choosing an exposure-outcome hypothesis to explore. e.g. LDL cholesterol on coronary heart disease. These data can be extracted from MR-Base:

```r
library(tryx)
a <- extract_instruments(300)
b <- extract_outcome_data(a$SNP, 7)
dat <- harmonise_data(a,b)
```

We can now perform the analysis:

```r
tryxscan <- tryx.scan(dat)
```

This will do the following:

1. Find outlier SNPs in the exposure-outcome analysis
2. Find traits in the MR-Base database that those outliers associate with. These traits are known as 'candidate traits'
3. Extract instruments for those 'candidate traits'
4. Perform MR of each of those 'candidate traits' against the exposure and the outcome

See the `?tryx.scan` for options on the parameters for this analysis. e.g. You can specify your own set of outliers, for example SNPs that have extreme p-values in the outcome GWAS

```r
x <- as.character(subset(dat, pval.outcome < 5e-8)$SNP)
tryxscan <- tryx.scan(dat, outliers=x)
```

The next steps are to determine which of the candidate traits are of interest (e.g. using p-value thresholds), visualise the results, and adjust the exposure-outcome estimates based on knowledge of the 'candidate trait' associations.

### Significant candidate traits

One can determine which of the putative associations might be 'interesting' in different ways. We have provided a simple convenience function to apply different multiple testing corrections. e.g.

```r
tryxscan <- tryx.sig(tryxscan)
```

Will by default use FDR of 5%. See `?tryx.sig` for more options.

### Visualisation

To produce a basic diagram of the connectivity of SNPs, candidate traits, exposure and outcome:

```r
tryx.network(tryxscan)
```

This shows that some candidate traits influence the exposure only, the outcome only, or both the exposure and the outcome. 

You can also create a volcano plot of the candidate-exposure and/or candidate-outcome associations. e.g. to show exposures and outcomes

```r
volcano_plot(
    rbind(tryxscan$candidate_exposure_mr, tryxscan$candidate_outcome_mr)
)
```

or e.g. just exposures

```r
volcano_plot(tryxscan$candidate_exposure_mr)
```

### Adjustment

Finally, to adjust the SNP effects on the exposure and outcome traits given their influences on the candidate traits, we can run:

```r
tryxanalysis <- tryx.analyse(tryxscan)
```

By default, this adjusts for the trait that has the largest impact for a particular SNP. 

The `tryxanalysis$estimates` show the adjusted effect estimates. A plot is generated showing how SNP effects have changed due to candidate trait adjustments in `tryxanalysis$plot`.

---


## To Do

### Plot

- categorise traits
- print SNP names or gene names
- option to have no names

### Scan

- Implement cooks distance as option for finding outliers
- Implement MR PRESSO as option for finding outliers

### Scan output

- calculate significant associations (currently in plot function - should remove from here)

### Analysis

- Simulate improvement in multiple testing correction when filtering by outlier associations
