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

---


## Implementing R6

1. specify data
2. find outliers
3. do all scans

OR

3. scan candidates
4. get instruments for candidates
5. get candidate instruments for outcome
6. get candidate instruments for exposure
7. get exposure instruments for candidates


8. harmonise all data
9. mr


GIB BROKE THIS

```
library(devtools)
install_github("mrcieu/TwoSampleMR@ieugwasr")
load_all()
a <- extract_instruments("ieu-a-300")
b <- extract_outcome_data(a$SNP, "ieu-a-7", access_token=NULL)
dat <- harmonise_data(a,b)

x <- Tryx$new()
x$input(dat)
x$get_outliers()
x$set_candidate_traits()

x$scan_candidate_traits() rename this to:
x$scan()

x$extractions(arg1,arg2,arg3,arg4) does these:
x$candidate_instruments(arg1)
x$outcome_instruments(arg2, arg3)
x$exposure_instruments(arg4)
x$exposure_candidate_instruments()

x$harmonise() does these:
x$candidate_outcome_dat()
x$candidate_exposure_dat()
x$exposure_candidate_dat()

x$perform_mr() rename this to:
x$mr()

# to do:
x$scan()

x$plots...
x$adjustments...
```


EXAMPLE OF CALLING A FUNCTION

```
Test <- R6::R6Class("test", list(
  output = 0
))

Test$set("public", "add", function(x) {
    self$output <- self$output + x
})


Test$set("public", "calladd", function(x) {
    self$add(x)
})

Test$set("public", "subt", function(x) self$output <- self$output - x)

test <- Test$new()

test$output
test$add(1)
test$calladd(1)
```


LAYOUT?

```
Tryx$dat        # input data
Tryx$output     # this is scan result
Tryx$output$outliers
           $candidate
           $mr_...


Tryx$adjustment
               $estimates
               $plot



```



## Guide

### Basic analysis

The following analyses should run within 10 minutes, depending on internet speed and the traffic that the MR-Base servers are experiencing.

Begin by choosing an exposure-outcome hypothesis to explore. e.g. LDL cholesterol on coronary heart disease. These data can be extracted from MR-Base:

```r
library(tryx)
a <- extract_instruments(300)
b <- extract_outcome_data(a$SNP, 7, access_token=NULL)
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


