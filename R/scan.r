

#' Outlier scan
#' 
#' A simple wrapper function.
#' Using a summary set, find outliers in the MR analysis between the pair of trais.
#' Find other 'candidate traits' associated with those outliers.
#' Perform MR of each of those candidate traits with the original exposure and outcome
#' 
#' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used
#' @param outliers Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed
#' @param search_threshold The p-value threshold for detecting an association between an outlier and a candidate trait. Default is 5e-8
#' @param id_list The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes()
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is strategy1. Can also provide basic MR methods e.g. mr_ivw, mr_weighted_mode etc
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
tryx.scan <- function(dat, outliers="RadialMR", use_proxies=FALSE, search_threshold=5e-8, id_list="default", include_outliers=FALSE, mr_method="strategy1")
{
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
		radial <- RadialMR::RadialMR(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP, "IVW", "YES", "NO", 0.05/nrow(dat), "NO")
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
			filter(! id %in% c(dat$id.exposure[1], dat$id.outcome[1]))
		id_list <- ids$id
		message("Using default list of ", nrow(ids), " traits")
	}
	output$id_list <- id_list

	output$search <- extract_outcome_data(radial$outliers$SNP, ids$id, proxies=use_proxies)
	output$search$sig <- output$search$pval.outcome < search_threshold
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

	message("Number of candidate - outcome associations: ", sum(tryxscan$candidate_outcome_mr$sig))
	message("Number of candidate - exposure associations: ", sum(tryxscan$candidate_exposure_mr$sig))
	message("Number of exposure - candidate associations: ", sum(tryxscan$exposure_candidate_mr$sig))
	return(tryxscan)
}


#' Analyse tryx results
#' 
#' This returns various heterogeneity statistics, IVW estimates for raw, 
#' adjusted and outlier removed datasets, and summary of peripheral 
#' traits detected etc.
#' 
#' @param tryxscan Output from \code{tryx.scan}
#' @param plot Whether to plot or not. Default is TRUE
#' @param filter_duplicate_outliers Whether to only allow each putative outlier to be adjusted by a single trait (in order of largest divergence). Default is TRUE.
tryx.analyse <- function(tryxscan, plot=TRUE, filter_duplicate_outliers=TRUE)
{

	analysis <- list()
	adj_full <- tryx.adjustment(tryxscan)
	if(nrow(adj_full) == 0)
	{
		return(NULL)
	}
	analysis$adj_full <- adj_full
	if(filter_duplicate_outliers)
	{
		adj <- adj_full %>% arrange(d) %>% filter(!duplicated(SNP))
	} else {
		adj <- adj_full
	}
	analysis$adj <- adj

	# Detection
	if("simulation" %in% names(tryxscan))
	{
		detection <- list()
		detection$nu1_correct <- sum(adj$SNP %in% 1:tryxscan$simulation$nu1 & adj$what != "p->y")
		detection$nu1_incorrect <- sum(adj$SNP %in% 1:tryxscan$simulation$nu1 & adj$what == "p->y")
		detection$nu2_correct <- sum(adj$SNP %in% (1:tryxscan$simulation$nu2 + tryxscan$simulation$nu1) & adj$what == "p->y")
		detection$nu2_incorrect <- sum(adj$SNP %in% (1:tryxscan$simulation$nu2 + tryxscan$simulation$nu1) & adj$what != "p->y")
		detection$no_outlier_flag <- tryxscan$simulation$no_outlier_flag
		analysis$detection <- detection
	}

	cpg <- require(ggrepel)
	if(!cpg)
	{
		stop("Please install the ggrepel package\ninstall.packages('ggrepel')")
	}

	dat <- subset(tryxscan$dat, mr_keep, select=c(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome))
	dat$ratio <- dat$beta.outcome / dat$beta.exposure
	dat$weights <- sqrt(dat$beta.exposure^2 / dat$se.outcome^2)
	dat$ratiow <- dat$ratio * dat$weights
	dat$what <- "Unadjusted"
	dat$candidate <- "NA"
	dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))

	analysis$Q$full_Q <- adj$Q[1]
	analysis$Q$num_reduced <- sum(adj$adj.qi < adj$qi)
	analysis$Q$mean_d <- mean(adj$d)


	temp <- subset(adj, select=c(SNP, adj.beta.exposure, adj.beta.outcome, adj.se.exposure, adj.se.outcome, candidate))
	temp$qi <- cochrans_q(temp$adj.beta.outcome / temp$adj.beta.exposure, temp$adj.se.outcome / abs(temp$adj.beta.exposure))
	ind <- grepl("\\|", temp$candidate)
	if(any(ind))
	{
		temp$candidate[ind] <- sapply(strsplit(temp$candidate[ind], split=" \\|"), function(x) x[[1]])
	}
	names(temp) <- gsub("adj.", "", names(temp))
	temp$what <- "Adjusted"
	temp$weights <- sqrt(temp$beta.exposure^2 / temp$se.outcome^2)
	temp$ratio <- temp$beta.outcome / temp$beta.exposure
	temp$ratiow <- temp$ratio * temp$weights

	dat_adj <- rbind(temp, dat) %>% filter(!duplicated(SNP))
	dat_adj$qi <- cochrans_q(dat_adj$beta.outcome / dat_adj$beta.exposure, dat_adj$se.outcome / abs(dat_adj$beta.exposure))
	analysis$Q$adj_Q <- sum(dat_adj$qi)

	analysis

	dat_rem <- subset(dat, !SNP %in% temp$SNP)

	est1 <- summary(lm(ratiow ~ -1 + weights, data=dat_rem))
	est2 <- summary(lm(ratiow ~ -1 + weights, data=dat))
	est3 <- summary(lm(ratiow ~ -1 + weights, data=dat_adj))

	estimates <- data.frame(
		est=c("Outliers removed", "Raw", "Outliers adjusted"),
		b=c(coefficients(est1)[1,1], coefficients(est2)[1,1], coefficients(est3)[1,1]), 
		se=c(coefficients(est1)[1,2], coefficients(est2)[1,2], coefficients(est3)[1,2]), 
		pval=c(coefficients(est1)[1,4], coefficients(est2)[1,4], coefficients(est3)[1,4]),
		nsnp=c(nrow(dat_rem), nrow(dat), nrow(dat_adj)),
		Q = c(sum(dat_rem$qi), sum(dat$qi), sum(dat_adj$qi)),
		int=0
	)
	estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q) 

	analysis$estimates <- estimates

	if(plot)
	{	
		temp2 <- merge(dat, temp, by="SNP")
		labs <- rbind(
			data.frame(label=temp2$SNP, x=temp2$weights.x, y=temp2$weights.x * temp2$ratio.x),
			data.frame(label=temp$candidate, x=temp$weights, y=temp$weights * temp$ratio)
		)
		p <- ggplot(rbind(dat, temp), aes(y=ratiow, x=weights)) +
		geom_abline(data=estimates, aes(slope=b, intercept=0, linetype=est)) +
		geom_label_repel(data=labs, aes(x=x, y=y, label=label), size=2, segment.color = "grey50") +
		geom_point(aes(colour=what)) +
		geom_segment(data=temp2, colour="grey50", aes(x=weights.x, xend=weights.y, y=ratiow.x, yend=ratiow.y), arrow = arrow(length = unit(0.01, "npc"))) +
		labs(colour="") +
		xlim(c(0, max(dat$weights))) +
		ylim(c(min(0, dat$ratiow, temp$ratiow), max(dat$ratiow, temp$ratiow)))
		analysis$plot <- p
	}
	return(analysis)
}


cochrans_q <- function(b, se)
{
	xw <- sum(b / se^2) / sum(1/se^2)
	qi <- (1/se^2) * (b - xw)^2
	return(qi)
}

bootstrap_path <- function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000)
{
	res <- rep(0, nboot)
	for(i in 1:nboot)
	{
		res[i] <- rnorm(1, gx, gx.se) - rnorm(1, gp, gp.se) * rnorm(1, px, px.se)
	}
	pe <- gx - gp * px
	return(c(pe, sd(res)))
}

radialmr <- function(dat, outlier=NULL)
{
	library(ggplot2)
	beta.exposure <- dat$beta.exposure
	beta.outcome <- dat$beta.outcome
	se.outcome <- dat$se.outcome
	w <- sqrt(beta.exposure^2 / se.outcome^2)
	ratio <- beta.outcome / beta.exposure
	ratiow <- ratio*w
	dat <- data.frame(w=w, ratio=ratio, ratiow=ratiow)
	if(is.null(outlier))
	{
		dat2 <- dat
	} else {
		dat2 <- dat[-c(outlier), ]
	}
	mod <- lm(ratiow ~ -1 + w, dat2)$coefficients[1]
	ggplot(dat, aes(x=w, y=ratiow)) +
	geom_point() +
	geom_abline(slope=mod, intercept=0)
}

