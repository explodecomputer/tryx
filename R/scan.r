#' Outlier scan
#' 
#' A simple wrapper function.
#' Using a summary set, find outliers in the MR analysis between the pair of trais.
#' Find other 'candidate traits' associated with those outliers.
#' Perform MR of each of those candidate traits with the original exposure and outcome
#' 
#' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.
#' @param outliers Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers.
#' @param outlier_correction Defualt = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#' @param outlier_threshold If outlier_correction = "none" then the p-value threshold for detecting outliers is by default 0.05.
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
#' @param search_correction Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#' @param search_threshold If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.
#' @param id_list The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes().
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is mr_ivw. Can also provide basic MR methods e.g. mr_weighted_mode, mr_weighted_median etc. Also possible to use "strategy1" which performs IVW in the first instance, but then weighted mode for associations with high heterogeneity.
#'
#' @export
#' @return List
#' dat  Cleaned dat input
#' radialmr  Results from RadialMR analysis
#' outliers  List of outliers used
#' id_list  List of GWAS IDs used
#' search  Result from search of outliers against GWAS IDs
#' candidate_instruments  Instruments for candidate traits
#' candidate_outcome  Extracted instrument SNPs from outcome
#' candidate_outcome_dat  Harmonised candidate - outcome dataset
#' candidate_outcome_mr  MR analysis of candidates against outcome
#' candidate_exposure   Extracted instrument SNPs from exposure
#' candidate_exposure_dat  Harmonised candidate - exposure dataset
#' candidate_exposure_mr  MR analysis of candidates against exposure
tryx.scan <- function(dat, outliers="RadialMR", outlier_correction="none", outlier_threshold=ifelse(outlier_correction=="none", 0.05/nrow(dat), 0.05), use_proxies=FALSE, search_correction="none", search_threshold=ifelse(search_correction=="none", 5e-8, 0.05), id_list="default", include_outliers=FALSE, mr_method="mr_ivw")
{

	stopifnot(search_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
	stopifnot(search_threshold > 0 & search_threshold < 1)

	stopifnot(outlier_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
	stopifnot(outlier_threshold > 0 & outlier_threshold < 1)

	# Get outliers

	output <- list()
	if(length(unique(dat$id.exposure)) > 1 | length(unique(dat$id.outcome)) > 1)
	{
		message("Warning! Multiple exposure/outcome combinations found")
		message("Only using first exposure / outcome combination")
	}
	dat <- subset(dat, id.exposure == id.exposure[1] & id.outcome == id.outcome[1] & mr_keep)
	output$dat <- dat

	stopifnot(length(mr_method) == 1)
	stopifnot(mr_method %in% mr_method_list()$obj | mr_method == "strategy1")

	if(outliers[1] == "RadialMR")
	{
		message("Using RadialMR package to detect outliers")
		cpg <- require(RadialMR)
		if(!cpg)
		{
			stop("Please install the RadialMR package\ndevtools::install_github('WSpiller/RadialMR')")
		}
		cpg <- require(dplyr)
		if(!cpg)
		{
			stop("Please install the RadialMR package\ndevtools::install_github('WSpiller/RadialMR')")
		}


		# radial <- RadialMR::ivw_radial(RadialMR::format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP), alpha=0.05/nrow(dat), weights=3)

		radialor <- RadialMR::ivw_radial(RadialMR::format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP), alpha=1, weights=3)

		# apply outlier_correction method with outlier_threshold to radial SNP-Q statistics
		radialor$data$rdpadj <- p.adjust(radialor$data$Qj_Chi, outlier_correction)
		radial <- radialor
		radial$outliers <- subset(radialor$data, radialor$data$rdpadj < outlier_threshold, select=c(SNP, Qj, Qj_Chi, rdpadj))
		colnames(radial$outliers) = c("SNP", "Q_statistic", "p.value", "adj.p.value")
		rownames(radial$outliers) <- 1:nrow(radial$outliers)

		# apply outlier_correction method with outlier_threshold to radial SNP-Q statistics


		if(radial$outliers[1] == "No significant outliers")
		{
			message("No outliers found")
			message("Try changing the outlier_threshold parameter")
			return(NULL)
		}
		outliers <- as.character(radial$outliers$SNP)
		message("Identified ", length(outliers), " outliers")
		output$radialmr <- radial
		if(length(outliers) == 0)
		{
			message("No outliers found. Exiting")
			return(output)
		}
	} else {
		nout <- length(outliers)
		outliers <- subset(dat, SNP %in% outliers)$SNP
		message(length(outliers), " of ", nout, " of specified outliers present in dat")
		if(length(outliers) == 0)
		{
			message("No outliers found. Exiting")
			return(output)
		}
	}
	output$outliers <- outliers

	# Find associations with outliers

	if(id_list[1] == "default")
	{
		ao <- suppressMessages(available_outcomes())
		ids <- subset(ao, priority == 1 & nsnp > 500000 & sample_size > 5000) %>%
			arrange(desc(sample_size)) %>%
			filter(!duplicated(trait), mr == 1) %>%
			filter(!author %in% c("Shin", "Kettunen", "Roederer")) %>%
			filter(!grepl("UKB-b", id)) %>%
			filter(! id %in% c(dat$id.exposure[1], dat$id.outcome[1]))
		id_list <- ids$id
		message("Using default list of ", nrow(ids), " traits")
	}
	output$id_list <- id_list

	output$search <- extract_outcome_data(outliers, output$id_list, proxies=use_proxies)
	padj <- p.adjust(output$search$pval.outcome, search_correction)
	output$search$sig <- padj < search_threshold
	out2 <- subset(output$search, sig)
	if(nrow(out2) == 0)
	{
		message("Outliers do not associate with any other traits. Try relaxing the search_threshold")
		return(output)
	}
	message("Found ", length(unique(out2$id.outcome)), " candidate traits associated with outliers at p < ", search_threshold)

	message("Finding instruments for candidate traits")
	output$candidate_instruments <- extract_instruments(unique(out2$id.outcome))

	if(nrow(output$candidate_instruments) == 0)
	{
		message("No instruments available for the candidate traits")
		return(output)
	}

	if(!include_outliers)
	{
		message("Removing outlier SNPs from candidate trait instrument lists")
		output$candidate_instruments <- group_by(output$candidate_instruments, id.exposure) %>%
			do({
				x <- .
				y <- subset(out2, id.outcome == x$id.exposure[1])
				x <- subset(x, !SNP %in% y$SNP)
				x
			})
	}

	if(nrow(output$candidate_instruments) == 0)
	{
		message("No instruments available for the candidate traits")
		return(output)
	}

	message(length(unique(output$candidate_instruments$id.exposure)), " traits with at least one instrument")

######

	message("Looking up candidate trait instruments for ", dat$outcome[1])
	output$candidate_outcome <- extract_outcome_data(unique(output$candidate_instruments$SNP), dat$id.outcome[1], proxies=use_proxies)
	if(is.null(output$candidate_outcome))
	{
		message("None of the candidate trait instruments available for ", dat$outcome[1])
		return(output)
	}
	message(nrow(output$candidate_outcome), " instruments extracted for ", dat$outcome[1])
	
	output$candidate_outcome_dat <- suppressMessages(harmonise_data(output$candidate_instruments, output$candidate_outcome))
	output$candidate_outcome_dat <- subset(output$candidate_outcome_dat, mr_keep)
	if(nrow(output$candidate_outcome_dat) == 0)
	{
		message("None of the candidate trait instruments available for ", dat$outcome[1], " after harmonising")
		return(output)
	}

	message("Performing MR of ", length(unique(output$candidate_outcome_dat$id.exposure)), " candidate traits against ", dat$outcome[1])
	if(mr_method == "strategy1")
	{
		temp <- strategy1(output$candidate_outcome_dat)
		output$candidate_outcome_mr <- temp$res
		output$candidate_outcome_mr_full <- temp
	} else {
		output$candidate_outcome_mr <- suppressMessages(mr(output$candidate_outcome_dat, method_list=c("mr_wald_ratio", mr_method)))
	}


######

	message("Looking up candidate trait instruments for ", dat$exposure[1])
	output$candidate_exposure <- extract_outcome_data(unique(output$candidate_instruments$SNP), dat$id.exposure[1], proxies=use_proxies)
	if(is.null(output$candidate_exposure))
	{
		message("None of the candidate trait instruments available for ", dat$exposure[1])
		return(output)
	}
	message(nrow(output$candidate_exposure), " instruments extracted for ", dat$exposure[1])
	
	output$candidate_exposure_dat <- suppressMessages(harmonise_data(output$candidate_instruments, output$candidate_exposure))
	output$candidate_exposure_dat <- subset(output$candidate_exposure_dat, mr_keep)
	if(nrow(output$candidate_exposure_dat) == 0)
	{
		message("None of the candidate trait instruments available for ", dat$exposure[1], " after harmonising")
		return(output)
	}

	message("Performing MR of ", length(unique(output$candidate_exposure_dat$id.exposure)), " candidate traits against ", dat$exposure[1])
	if(mr_method == "strategy1")
	{
		temp <- strategy1(output$candidate_exposure_dat)
		output$candidate_exposure_mr <- temp$res
		output$candidate_exposure_mr_full <- temp
	} else {
		output$candidate_exposure_mr <- suppressMessages(mr(output$candidate_exposure_dat, method_list=c("mr_wald_ratio", mr_method)))
	}

#####

	message("Looking up exposure instruments for ", length(unique(out2$id.outcome)), " candidate traits")
	output$exposure_candidate <- extract_outcome_data(unique(output$dat$SNP), unique(out2$id.outcome), proxies=use_proxies)
	if(is.null(output$exposure_candidate))
	{
		message("None of the candidate trait instruments available for ", dat$exposure[1])
		return(output)
	}
	message(nrow(output$candidate_exposure), " instruments extracted")
	
	temp <- subset(output$dat, select=c(SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, id.exposure, exposure, units.exposure))

	if(!include_outliers)
	{
		message("Removing outlier SNPs from candidate trait outcome lists")
		n1 <- nrow(output$exposure_candidate)
		output$exposure_candidate <- group_by(output$exposure_candidate, id.outcome) %>%
			do({
				x <- .
				y <- subset(out2, id.outcome == x$id.outcome[1])
				x <- subset(x, !SNP %in% y$SNP)
				x
			})
		message("Removed ", n1 - nrow(output$exposure_candidate), " outlier SNPs")
	}


	output$exposure_candidate_dat <- suppressMessages(harmonise_data(temp, output$exposure_candidate))
	output$exposure_candidate_dat <- subset(output$exposure_candidate_dat, mr_keep)
	if(nrow(output$exposure_candidate_dat) == 0)
	{
		message("None of the candidate traits have the ", dat$exposure[1], " instruments after harmonising")
		return(output)
	}

	message("Performing MR of ", dat$exposure[1], " against ", length(unique(output$exposure_candidate_dat$id.outcome)), " candidate traits")
	if(mr_method == "strategy1")
	{
		temp <- strategy1(output$exposure_candidate_dat)
		output$exposure_candidate_mr <- temp$res
		output$exposure_candidate_mr_full <- temp
	} else {
		output$exposure_candidate_mr <- suppressMessages(mr(output$exposure_candidate_dat, method_list=c("mr_wald_ratio", mr_method)))
	}

	return(output)
}


#' MR Strategy 1
#' 
#' How to choose the result for a set of different MR analysies?
#' Simple strategy:
#' Use Wald ratio if only one SNP
#' Use IVW if more than one SNP and heterogeneity is low
#' Use weighted mode if more than some minimum number of SNPs and heterogeneity is high
#' 
#' @param dat Output from harmonise_data function
#' @param het_threshold The p-value threshold for Cochran's Q - if lower than this threshold then run weighted mode. Default p = 0.05
#' @param ivw_max_snp Maximum SNPs to allow IVW result even if heterogeneity is high. Default = 1
strategy1 <- function(dat, het_threshold=0.05, ivw_max_snp=1)
{
	message("First pass: running ", length(unique(paste(dat$id.exposure, dat$id.outcome))), " analyses with Wald ratio or IVW")
	a <- suppressMessages(mr(dat, method_list=c("mr_ivw", "mr_wald_ratio")))
	b <- subset(a, method == "Wald ratio")
	message(nrow(b), " analyses can only run Wald ratio")
	c <- subset(a, method == "Inverse variance weighted")
	message(nrow(c), " analyses can run IVW or mode")
	if(nrow(c) > 0)
	{
		d <- suppressMessages(mr_heterogeneity(dat, method_list="mr_ivw"))
		rerun <- subset(d, Q_pval < het_threshold & Q_df >= (ivw_max_snp-1))
		rerun <- paste(rerun$id.exposure, rerun$id.outcome)
		if(length(rerun) < 0)
		{
			message("All eligible IVW results have low heterogeneity")
			return(list(res=a, all=NULL, heterogeneity=d))
		}
		e <- subset(c, !paste(id.exposure, id.outcome) %in% rerun)
		message("Rerunning MR with weighted modal estimator for ", length(rerun), " analysis - this may be quite slow")
		f <- suppressMessages(mr(
			subset(dat, paste(id.exposure, id.outcome) %in% rerun),
			method_list="mr_weighted_mode"
		))
		res <- rbind(b, e, f)
		return(list(res=res, all=a, heterogeneity=d))
	} else {
		return(list(res=a, all=NULL, heterogeneity=NULL))
	}
}


#' Identify putatively significant associations in the outlier scan
#' 
#' @param mr_threshold_method This is the argument to be passed to \code{p.adjust}. Default is "fdr". If no p-value adjustment is to be applied then specify "unadjusted"
#' @param mr_threshold Threshold to declare significance
#' @export
#' @return Same as outlier_scan but the candidate_exposure_mr and candidate_outcome_mr objects have an extra pval_adj and sig column each
tryx.sig <- function(tryxscan, mr_threshold_method = "fdr", mr_threshold = 0.05)
{
	stopifnot("candidate_outcome_mr" %in% names(tryxscan))
	stopifnot("candidate_exposure_mr" %in% names(tryxscan))

	# Use threshold to retain causal relationships
	stopifnot(length(mr_threshold_method) == 1)
	if(mr_threshold_method != "unadjusted")
	{
		message("Adjusting p-value")
		tryxscan$candidate_outcome_mr$pval_adj <- p.adjust(tryxscan$candidate_outcome_mr$pval, mr_threshold_method)
		tryxscan$candidate_outcome_mr$sig <- tryxscan$candidate_outcome_mr$pval_adj < mr_threshold
		tryxscan$candidate_exposure_mr$pval_adj <- p.adjust(tryxscan$candidate_exposure_mr$pval, mr_threshold_method)
		tryxscan$candidate_exposure_mr$sig <- tryxscan$candidate_exposure_mr$pval_adj < mr_threshold
		tryxscan$exposure_candidate_mr$pval_adj <- p.adjust(tryxscan$exposure_candidate_mr$pval, mr_threshold_method)
		tryxscan$exposure_candidate_mr$sig <- tryxscan$exposure_candidate_mr$pval_adj < mr_threshold
	} else {
		tryxscan$candidate_outcome_mr$sig <- tryxscan$candidate_outcome_mr$pval < mr_threshold
		tryxscan$candidate_exposure_mr$sig <- tryxscan$candidate_exposure_mr$pval < mr_threshold
		tryxscan$exposure_candidate_mr$sig <- tryxscan$exposure_candidate_mr$pval < mr_threshold
	}

	message("* * * *")
	message("Number of candidate - outcome associations: ", sum(tryxscan$candidate_outcome_mr$sig))
	message("* * * *")
	message(paste(subset(tryxscan$candidate_outcome_mr, sig)$exposure, collapse="\n"))
	message("")
	message("* * * *")
	message("Number of candidate - exposure associations: ", sum(tryxscan$candidate_exposure_mr$sig))
	message("* * * *")
	message(paste(subset(tryxscan$candidate_exposure_mr, sig)$exposure, collapse="\n"))
	message("")
	message("* * * *")
	message("Number of exposure - candidate associations: ", sum(tryxscan$exposure_candidate_mr$sig))
	message("* * * *")
	message(paste(subset(tryxscan$candidate_exposure_mr, sig)$exposure, collapse="\n"))
	return(tryxscan)
}


