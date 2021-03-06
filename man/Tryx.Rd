% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tryxclass.R
\name{Tryx}
\alias{Tryx}
\title{Class for MR-TRYX analysis}
\description{
A simple wrapper function.
Using a summary set, find outliers in the MR analysis between the pair of trais.
Find other 'candidate traits' associated with those outliers.
Perform MR of each of those candidate traits with the original exposure and outcome.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Tryx$new()}}
\item \href{#method-print}{\code{Tryx$print()}}
\item \href{#method-get_outliers}{\code{Tryx$get_outliers()}}
\item \href{#method-set_candidate_traits}{\code{Tryx$set_candidate_traits()}}
\item \href{#method-scan}{\code{Tryx$scan()}}
\item \href{#method-candidate_instruments}{\code{Tryx$candidate_instruments()}}
\item \href{#method-outcome_instruments}{\code{Tryx$outcome_instruments()}}
\item \href{#method-exposure_instruments}{\code{Tryx$exposure_instruments()}}
\item \href{#method-exposure_candidate_instruments}{\code{Tryx$exposure_candidate_instruments()}}
\item \href{#method-extractions}{\code{Tryx$extractions()}}
\item \href{#method-candidate_outcome_dat}{\code{Tryx$candidate_outcome_dat()}}
\item \href{#method-candidate_exposure_dat}{\code{Tryx$candidate_exposure_dat()}}
\item \href{#method-exposure_candidate_dat}{\code{Tryx$exposure_candidate_dat()}}
\item \href{#method-harmonise}{\code{Tryx$harmonise()}}
\item \href{#method-mr}{\code{Tryx$mr()}}
\item \href{#method-mrtryx}{\code{Tryx$mrtryx()}}
\item \href{#method-tryx.sig}{\code{Tryx$tryx.sig()}}
\item \href{#method-adjustment}{\code{Tryx$adjustment()}}
\item \href{#method-adjustment.mv}{\code{Tryx$adjustment.mv()}}
\item \href{#method-analyse}{\code{Tryx$analyse()}}
\item \href{#method-analyse.mv}{\code{Tryx$analyse.mv()}}
\item \href{#method-manhattan_plot}{\code{Tryx$manhattan_plot()}}
\item \href{#method-clone}{\code{Tryx$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\subsection{Method \code{new()}}{
Create a new dataset and initialise an R interface
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$new(dat)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Dataset from TwoSampleMR::harmonise_data}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$print(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_outliers"></a>}}
\subsection{Method \code{get_outliers()}}{
Detect outliers in exposure-outcome dataset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$get_outliers(
  dat = self$output$dat,
  outliers = "RadialMR",
  outlier_correction = "none",
  outlier_threshold = ifelse(outlier_correction == "none", 0.05/nrow(dat), 0.05)
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.}

\item{\code{outliers}}{Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers.}

\item{\code{outlier_correction}}{Defualt = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").}

\item{\code{outlier_threshold}}{If outlier_correction = "none" then the p-value threshold for detecting outliers is by default 0.05.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_candidate_traits"></a>}}
\subsection{Method \code{set_candidate_traits()}}{
Set a list of GWAS IDs used.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$set_candidate_traits(id_list = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id_list}}{The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes().}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scan"></a>}}
\subsection{Method \code{scan()}}{
Search for candidate traits associated with outliers.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$scan(
  dat = self$output$dat,
  search_correction = "none",
  search_threshold = ifelse(search_correction == "none", 5e-08, 0.05),
  use_proxies = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}

\item{\code{search_correction}}{Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none.}

\item{\code{search_threshold}}{If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.}

\item{\code{use_proxies}}{Whether to use proxies when looking up associations. FALSE by default for speed.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-candidate_instruments"></a>}}
\subsection{Method \code{candidate_instruments()}}{
Obtain instruments for the candidate traits.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$candidate_instruments(
  candidate_instruments = NULL,
  include_outliers = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{candidate_instruments}}{Instruments for candidate traits.}

\item{\code{include_outliers}}{When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-outcome_instruments"></a>}}
\subsection{Method \code{outcome_instruments()}}{
Extract instrument for candidate trait instruments for the original outcome.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$outcome_instruments(
  candidate_outcome = NULL,
  dat = self$output$dat,
  use_proxies = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{candidate_outcome}}{Extracted instrument SNPs from outcome.}

\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}

\item{\code{use_proxies}}{Whether to use proxies when looking up associations. FALSE by default for speed.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exposure_instruments"></a>}}
\subsection{Method \code{exposure_instruments()}}{
Extract instrument for candidate trait instruments for the original exposure.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$exposure_instruments(
  candidate_exposure = NULL,
  dat = self$output$dat,
  use_proxies = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{candidate_exposure}}{Extracted instrument SNPs from exposure.}

\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}

\item{\code{use_proxies}}{Whether to use proxies when looking up associations. FALSE by default for speed.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exposure_candidate_instruments"></a>}}
\subsection{Method \code{exposure_candidate_instruments()}}{
Extract instrument for the original exposure for the candidate traits.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$exposure_candidate_instruments(
  exposure_candidate = NULL,
  dat = self$output$dat,
  use_proxies = FALSE,
  include_outliers = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{exposure_candidate}}{Extracted instrument SNPs from exposure.}

\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}

\item{\code{use_proxies}}{Whether to use proxies when looking up associations. FALSE by default for speed.}

\item{\code{include_outliers}}{When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-extractions"></a>}}
\subsection{Method \code{extractions()}}{
Extract instruments for MR analyses.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$extractions(
  dat = self$output$dat,
  candidate_instruments = NULL,
  candidate_outcome = NULL,
  candidate_exposure = NULL,
  exposure_candidate = NULL,
  include_outliers = FALSE,
  use_proxies = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}

\item{\code{candidate_instruments}}{Instruments for candidate traits.}

\item{\code{candidate_outcome}}{Extracted instrument SNPs from outcome.}

\item{\code{candidate_exposure}}{Extracted instrument SNPs from exposure.}

\item{\code{exposure_candidate}}{Extracted instrument SNPs from exposure.}

\item{\code{include_outliers}}{When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.}

\item{\code{use_proxies}}{Whether to use proxies when looking up associations. FALSE by default for speed.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-candidate_outcome_dat"></a>}}
\subsection{Method \code{candidate_outcome_dat()}}{
Make a dataset for the candidate traits and the original outcome.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$candidate_outcome_dat(dat = self$output$dat)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-candidate_exposure_dat"></a>}}
\subsection{Method \code{candidate_exposure_dat()}}{
Make a dataset for the candidate traits and the original exposure.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$candidate_exposure_dat(dat = self$output$dat)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exposure_candidate_dat"></a>}}
\subsection{Method \code{exposure_candidate_dat()}}{
Make a dataset for the original exposure and the candidate traits.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$exposure_candidate_dat(dat = self$output$dat)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-harmonise"></a>}}
\subsection{Method \code{harmonise()}}{
Harmonised exposure - outcome dataset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$harmonise(dat = self$output$dat)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mr"></a>}}
\subsection{Method \code{mr()}}{
Perform MR anlayses of 1) candidate traits-outcome 2) candidate traits-exposure 3) exposure-candidate traits.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$mr(dat = self$output$dat, mr_method = "mr_ivw")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from TwoSampleMR::harmonise_data.}

\item{\code{mr_method}}{Method to use for candidate trait - exposure/outcome analysis. Default is mr_ivw. Can also provide basic MR methods e.g. mr_weighted_mode, mr_weighted_median etc. Also possible to use "strategy1" which performs IVW in the first instance, but then weighted mode for associations with high heterogeneity.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mrtryx"></a>}}
\subsection{Method \code{mrtryx()}}{
All-in-one: 1) Detect outlier 2) Search candidate traits 3) Perform MR of candidate traits and the outcome / exposure.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$mrtryx(
  dat = self$output$dat,
  outliers = "RadialMR",
  outlier_correction = "none",
  outlier_threshold = ifelse(outlier_correction == "none", 0.05/nrow(dat), 0.05),
  use_proxies = FALSE,
  search_correction = "none",
  search_threshold = ifelse(search_correction == "none", 5e-08, 0.05),
  include_outliers = FALSE,
  mr_method = "mr_ivw"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dat}}{Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.}

\item{\code{outliers}}{Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers.}

\item{\code{outlier_correction}}{Defualt = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").}

\item{\code{outlier_threshold}}{If outlier_correction = "none" then the p-value threshold for detecting outliers is by default 0.05.}

\item{\code{use_proxies}}{Whether to use proxies when looking up associations. FALSE by default for speed.}

\item{\code{search_correction}}{Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").}

\item{\code{search_threshold}}{If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.}

\item{\code{include_outliers}}{When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.}

\item{\code{mr_method}}{Method to use for candidate trait - exposure/outcome analysis. Default is mr_ivw. Can also provide basic MR methods e.g. mr_weighted_mode, mr_weighted_median etc. Also possible to use "strategy1" which performs IVW in the first instance, but then weighted mode for associations with high heterogeneity.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-tryx.sig"></a>}}
\subsection{Method \code{tryx.sig()}}{
Identify putatively significant associations in the outlier scan.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$tryx.sig(mr_threshold_method = "fdr", mr_threshold = 0.05)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mr_threshold_method}}{This is the argument to be passed to \code{p.adjust}. Default is "fdr". If no p-value adjustment is to be applied then specify "unadjusted".}

\item{\code{mr_threshold}}{Threshold to declare significance}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-adjustment"></a>}}
\subsection{Method \code{adjustment()}}{
Outlier adjustment estimation - How much of the heterogeneity due to the outlier can be explained by alternative pathways?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$adjustment(tryxscan = self$output, id_remove = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{tryxscan}}{Output from \code{x$mrtryx()}}

\item{\code{id_remove}}{List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.}

\item{\code{dat}}{Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-adjustment.mv"></a>}}
\subsection{Method \code{adjustment.mv()}}{
Similar to adjusment, but when there are multiple traits associated with a single variant.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$adjustment.mv(
  tryxscan = self$output,
  lasso = TRUE,
  id_remove = NULL,
  proxies = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{tryxscan}}{Output from \code{x$scan()}}

\item{\code{lasso}}{Whether to shrink the estimates of each trait within SNP. Default=TRUE.}

\item{\code{id_remove}}{List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.}

\item{\code{proxies}}{Look for proxies in the MVMR methods. Default = FALSE.}

\item{\code{dat}}{Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-analyse"></a>}}
\subsection{Method \code{analyse()}}{
This returns various heterogeneity statistics, IVW estimates for raw, adjusted and outlier removed datasets, and summary of peripheral traits detected etc.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$analyse(
  tryxscan = self$output,
  plot = TRUE,
  id_remove = NULL,
  filter_duplicate_outliers = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{tryxscan}}{Output from \code{x$scan()}.}

\item{\code{plot}}{Whether to plot or not. Default is TRUE.}

\item{\code{id_remove}}{List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.}

\item{\code{duplicate_outliers_method}}{Sometimes more than one trait will associate with a particular outlier. TRUE = only keep the trait that has the biggest influence on heterogeneity.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-analyse.mv"></a>}}
\subsection{Method \code{analyse.mv()}}{
Similar to analyse, but when there are multiple traits associated with a single variant.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$analyse.mv(
  tryxscan = self$output,
  lasso = TRUE,
  plot = TRUE,
  id_remove = NULL,
  proxies = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{tryxscan}}{Output from \code{x$scan()}}

\item{\code{lasso}}{Whether to shrink the estimates of each trait within SNP. Default=TRUE.}

\item{\code{id_remove}}{List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.}

\item{\code{proxies}}{Look for proxies in the MVMR methods. Default = FALSE.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-manhattan_plot"></a>}}
\subsection{Method \code{manhattan_plot()}}{
Draw a Manhattan style plot for candidate traits-outcome/exposure associations.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$manhattan_plot(
  what = "outcome",
  id_remove = NULL,
  y_scale = NULL,
  label = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{what}}{Analyse candidate-exposure ('exposure') or candidate-outcome ('outcome') associations}

\item{\code{id_remove}}{List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.}

\item{\code{y_scale}}{The scaling function to be applied to y scale.}

\item{\code{label}}{Display the names of the traits on the graph.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tryx$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
