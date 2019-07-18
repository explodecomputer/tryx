#' Outlier adjustment estimation
#' 
#' How much of the heterogeneity due to the outlier can be explained by alternative pathways?
#' 
#' @param tryxscan Output from \code{tryx.scan}
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#' @param duplicate_outliers_method Sometimes more than one trait will associate with a particular outlier. The options are as follows. heterogeneity: For every SNP, adjust for each potential p independently and keep just the best one based on the biggest reduction in heterogeneity. none: For every SNP, adjust for each potential p independently and keep  Now keep all of them. mv: For each SNP, adjust for all potential p jointly by performing multivariable MR to obtain p->y causal estimates. mv_lasso: For each SNP, adjust for all potential p jointly by selecting p based on lasso, and then using multivariable MR on remaining variables to obtain p->y causal estimates. lasso_selection: For each SNP, adjust for all potential p jointly by selecting p based on lasso and using the original estimates. lasso: For each SNP, adjust for all potential p jointly by selecting p based on lasso and using the lasso p->y estimates with the original standard errors. Default is 'heterogeneity'.

#'
#' @export
#' @return data frame of adjusted effect estimates and heterogeneity stats
tryx.adjustment <- function(tryxscan, id_remove=NULL, duplicate_outliers_method=c("heterogeneity", "none", "mv", "mv_lasso", "lasso_selection", "lasso")[1])
{

	if(!any(tryxscan$search$sig))
	{
		return(tibble())
	}


	# Get all CT related to outliers that are significant, and remove any IDs to be removed
	sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)

	# These are exposure MR results that are significant with removed IDs
	sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)

	# These are exposure MR results that are significant with removed IDs
	sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)

	# These original exposure outcome data
	dat <- tryxscan$dat

	# Add heterogeneity
	dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
	dat$Q <- sum(dat$qi)


	# perform adjustment within every SNP

	search_result <- subset(sig, id.outcome %in% sigo$id.exposure) %>%
		group_by(SNP) %>% 
		mutate(snpcount=n()) %>% 
		arrange(desc(snpcount), SNP)

	snplist <- unique(search_result$SNP)

	l <- list()
	if(duplicate_outliers_method %in% c("none", "heterogeneity"))
	{
		for(i in 1:nrow(sig))
		{
			l[[i]] <- adjust_each_entry(i, dat, sig, sigo)
		}

	} else {

		mvo <- perform_mv(search_result, snplist, sig, dat$id.exposure[1], dat$id.outcome[1])
		return(mvo)

		if(duplicate_outliers_method == "mv")
		{
			for(i in 1:length(mvo))
			{
				l[[i]] <- adjust_each_entry(i, dat, sig, mvo$mvo[[i]]$sigo_mv)
			}
		} else if(duplicate_outliers_method == "mv_lasso") {

			for(i in 1:length(mvo))
			{
				l[[i]] <- adjust_each_entry(i, dat, mvo$sig[[i]], mvo$mvo[[i]]$sigo_selected)
			}
			adjusted <- bind_rows(l)
			adjusted$d <- (adjusted$adj.qi / adjusted$adj.Q) / (adjusted$qi / adjusted$Q)

		} else if(duplicate_outliers_method == "lasso_selection") {

			for(i in 1:length(mvo))
			{
				l[[i]] <- adjust_each_entry(i, dat, mvo[[i]]$sig, sig)
			}

		} else if(duplicate_outliers_method == "lasso") {

			for(i in 1:length(mvo))
			{
				l[[i]] <- adjust_each_entry(i, dat, mvo[[i]]$sig, mvo[[i]]$mvo$sigo_lasso)
			}
		}
	}

	adjusted <- bind_rows(l)
	adjusted$d <- (adjusted$adj.qi / adjusted$adj.Q) / (adjusted$qi / adjusted$Q)

	if(duplicate_outliers_method == "heterogeneity")
	{
		adjusted <- adjusted %>% arrange(d) %>% filter(!duplicated(SNP))
	}
	return(adjusted)
}



perform_mv <- function(search_result, snplist, sig, id.exposure, id.outcome)
{
	mvo <- list()
	ssig <- list()
	for(i in 1:length(snplist))
	{
		mvo[[i]] <- list()
		message("Estimating joint effects of the following traits associated with ", snplist[i])
		temp <- subset(search_result, SNP %in% snplist[i])
		message(paste(unique(temp$outcome), collapse="\n"))
		message("Extracting instruments")
		mvexp <- suppressMessages(mv_extract_exposures(c(id.exposure, temp$id.outcome), find_proxies=proxies))
		message("Extracting outcome data for new instruments")
		mvout <- suppressMessages(extract_outcome_data(mvexp$SNP, id.outcome))
		mvdat <- suppressMessages(mv_harmonise_data(mvexp, mvout))
		message("Performing LASSO")
		b <- glmnet::cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
		c <- coef(b, s = "lambda.min")
		keeplist <- unique(c(rownames(c)[!c[,1] == 0], tryxscan$dat$exposure[1]))
		message("After shrinkage keeping:")
		message(paste(keeplist, collapse="\n"))
		mvexp2 <- subset(mvexp, exposure %in% keeplist)
		remsnp <- group_by(mvexp2, SNP) %>% summarise(mp = min(pval.exposure)) %>% filter(mp > 5e-8) %$% as.character(SNP)
		mvexp2 <- subset(mvexp2, !SNP %in% remsnp)
		mvout2 <- subset(mvout, SNP %in% mvexp2$SNP)
		mvdat2 <- suppressMessages(mv_harmonise_data(mvexp2, mvout2))

		mvo[[i]]$mv_multiple <- mv_multiple(mvdat)$result
		mvo[[i]]$mv_selection <- mv_multiple(mvdat2)$result
		mvo[[i]]$lasso <- coef(b) %>% as.matrix %>% {data_frame(exposure=rownames(.), b=.[,1])} %>% filter(exposure != "(Intercept)")
		remlist <- unique(mvo[[i]]$lasso$exposure[mvo[[i]]$lasso$b == 0])

		# Filter out sig to keep only selected traits per SNP
		ssig[[i]] <- subset(sig, !(SNP == snplist[i] & outcome %in% remlist))

		# Update sigo with lasso mv betas
		sigo_selected <- subset(sigo, exposure %in% mvo[[i]]$mv_selection$exposure)
		index <- match(sigo_selected$exposure, mvo[[i]]$mv_selection$exposure)
		stopifnot(all(sigo_selected$exposure == mvo[[i]]$mv_selection$exposure[index]))
		sigo_selected$b <- mvo[[i]]$mv_selection$b[index]
		sigo_selected$se <- mvo[[i]]$mv_selection$se[index]

		# Update sigo with mv betas
		sigo_mv <- subset(sigo, exposure %in% mvo[[i]]$mv_multiple$exposure)
		index <- match(sigo_mv$exposure, mvo[[i]]$mv_multiple$exposure)
		stopifnot(all(sigo_mv$exposure == mvo[[i]]$mv_multiple$exposure[index]))
		sigo_mv$b <- mvo[[i]]$mv_multiple$b[index]
		sigo_mv$se <- mvo[[i]]$mv_multiple$se[index]

		# Update sigo with lasso betas
		sigo_lasso <- subset(sigo, exposure %in% mvo[[i]]$lasso$exposure)
		index <- match(sigo_lasso$exposure, mvo[[i]]$lasso$exposure)
		stopifnot(all(sigo_lasso$exposure == mvo[[i]]$lasso$exposure[index]))
		sigo_lasso$b <- mvo[[i]]$lasso$b[index]

		# Gather
		mvo[[i]]$sigo_selected <- sigo_selected
		mvo[[i]]$sigo_mv <- sigo_mv
		mvo[[i]]$sigo_lasso <- sigo_lasso

	}
	return(list(sig=ssig, mvo=mvo))
}


adjust_each_entry <- function(i, dat, sig, sigo)
{
	a <- subset(dat, SNP == sig$SNP[i], select=c(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome, qi, Q))
	if(sig$id.outcome[i] %in% sigo$id.exposure)
	{
		a$candidate <- sig$outcome[i]
		a$i <- i
		message("   p->y:  ", a$SNP, "\t- ", sig$outcome[i])
		a$what <- "p->y"
		a$candidate.beta.exposure <- NA
		a$candidate.se.exposure <- NA
		a$candidate.beta.outcome <- sigo$b[sigo$id.exposure == sig$id.outcome[i]]
		a$candidate.se.outcome <- sigo$se[sigo$id.exposure == sig$id.outcome[i]]
		b <- bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
		a$adj.beta.exposure <- a$beta.exposure
		a$adj.se.exposure <- a$se.exposure
		a$adj.beta.outcome <- b[1]
		a$adj.se.outcome <- b[2]
		temp <- dat
		temp$beta.exposure[temp$SNP == a$SNP] <- a$adj.beta.exposure
		temp$se.exposure[temp$SNP == a$SNP] <- a$adj.se.exposure
		temp$beta.outcome[temp$SNP == a$SNP] <- a$adj.beta.outcome
		temp$se.outcome[temp$SNP == a$SNP] <- a$adj.se.outcome
		temp$qi <- cochrans_q(temp$beta.outcome / temp$beta.exposure, temp$se.outcome / abs(temp$beta.exposure))
		a$adj.qi <- temp$qi[temp$SNP == a$SNP]
		a$adj.Q <- sum(temp$qi)

		return(a)
	} else {
		return(NULL)
	}
}


#' Analyse tryx results
#' 
#' This returns various heterogeneity statistics, IVW estimates for raw, 
#' adjusted and outlier removed datasets, and summary of peripheral 
#' traits detected etc.
#' 
#' @param tryxscan Output from \code{tryx.scan}
#' @param plot Whether to plot or not. Default is TRUE
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#' @param filter_duplicate_outliers Whether to only allow each putative outlier to be adjusted by a single trait (in order of largest divergence). Default is TRUE.
#' 
#' @export
#' @return List of 
#' - adj_full: data frame of SNP adjustments for all candidate traits
#' - adj: The results from adj_full selected to adjust the exposure-outcome model
#' - Q: Heterogeneity stats
#' - estimates: Adjusted and unadjested exposure-outcome effects
#' - plot: Radial plot showing the comparison of different methods and the changes in SNP effects ater adjustment
tryx.analyse <- function(tryxscan, plot=TRUE, id_remove=NULL, duplicate_outliers_method=c("heterogeneity", "none", "mv", "mv_lasso", "lasso_selection", ""))
{

	cpg <- require(ggrepel)
	if(!cpg)
	{
		stop("Please install the ggrepel package\ninstall.packages('ggrepel')")
	}

	analysis <- list()
	adj <- tryx.adjustment(tryxscan, id_remove, duplicate_outliers_method)
	if(nrow(adj) == 0)
	{
		return(NULL)
	}
	analysis$adj <- adj

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

	dat_rem <- subset(dat, !SNP %in% temp$SNP)
	dat_rem2 <- subset(dat, !SNP %in% tryxscan$outliers)
	# adjust everything that's possible remove remaining outliers
	dat_rem3 <- rbind(temp, dat_rem2) %>% filter(!duplicated(SNP))

	est1 <- summary(lm(ratiow ~ -1 + weights, data=dat_rem))
	est2 <- summary(lm(ratiow ~ -1 + weights, data=dat))
	est3 <- summary(lm(ratiow ~ -1 + weights, data=dat_adj))
	est4 <- try(summary(lm(ratiow ~ -1 + weights, data=dat_rem2)))
	est5 <- try(summary(lm(ratiow ~ -1 + weights, data=dat_rem3)))
	if(class(est4) == "try-error")
	{
		est4 <- list(coefficients = matrix(NA, 2,4))
	}

	estimates <- data.frame(
		est=c("Raw", "Outliers removed (candidates)", "Outliers removed (all)", "Outliers adjusted", "Mixed", "Multivariable MR"),
		b=c(coefficients(est2)[1,1], coefficients(est1)[1,1], coefficients(est4)[1,1], coefficients(est3)[1,1], coefficients(est5)[1,1], tryxscan$mvres$result$b[1]), 
		se=c(coefficients(est2)[1,2], coefficients(est1)[1,2], coefficients(est4)[1,2], coefficients(est3)[1,2], coefficients(est5)[1,2], tryxscan$mvres$result$se[1]), 
		pval=c(coefficients(est2)[1,4], coefficients(est1)[1,4], coefficients(est4)[1,4], coefficients(est3)[1,4], coefficients(est5)[1,4], tryxscan$mvres$result$pval[1]),
		nsnp=c(nrow(dat), nrow(dat_rem), nrow(dat_rem2), nrow(dat_adj), nrow(dat_rem3), tryxscan$mvres$result$nsnp[1]),
		Q = c(sum(dat$qi), sum(dat_rem$qi), sum(dat_rem2$qi), sum(dat_adj$qi), sum(dat_rem3$qi), NA),
		int=0
	)
	estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q) 

	if("true_outliers" %in% names(tryxscan))
	{
		dat_true <- subset(dat, !SNP %in% tryxscan$true_outliers)
		est6 <- try(summary(lm(ratiow ~ -1 + weights, data=dat_true)))
		estimates <- bind_rows(estimates,
			data.frame(
				est=c("Oracle"),
				b=c(coefficients(est6)[1,1]), 
				se=c(coefficients(est6)[1,2]), 
				pval=c(coefficients(est6)[1,4]),
				nsnp=c(nrow(dat_true)),
				Q = c(sum(dat_true$qi)),
				int=0
			)
		)
	}
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
		geom_abline(data=estimates, aes(slope=b, intercept=0, colour=est)) +
		geom_label_repel(data=labs, aes(x=x, y=y, label=label), size=2, segment.color = "grey10") +
		geom_point() +
		geom_segment(data=temp2, colour="grey50", aes(x=weights.x, xend=weights.y, y=ratiow.x, yend=ratiow.y), arrow = arrow(length = unit(0.02, "npc"))) +
		labs(colour="") +
		xlim(c(0, max(dat$weights))) +
		ylim(c(min(0, dat$ratiow, temp$ratiow), max(dat$ratiow, temp$ratiow))) +
		scale_colour_brewer(type="qual")
		analysis$plot <- p
	}
	return(analysis)
}



#' Cochran's Q statistic
#' 
#' @param b vector of effecti 
#' @param se vector of standard errors
#' 
#' @return q values
#' 
#' @export
cochrans_q <- function(b, se)
{
	xw <- sum(b / se^2) / sum(1/se^2)
	qi <- (1/se^2) * (b - xw)^2
	return(qi)
}

bootstrap_path1 <- function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000)
{
	res <- rnorm(nboot, gx, gx.se) - rnorm(nboot, gp, gp.se) * rnorm(nboot, px, px.se)
	pe <- gx - gp * px
	return(c(pe, sd(res)))
}

bootstrap_path <- function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000)
{
	nalt <- length(gp)
	altpath <- data_frame(
		p = rnorm(nboot * nalt, gp, gp.se) * rnorm(nboot * nalt, px, px.se),
		b = rep(1:nboot, each=nalt)
	)
	altpath <- group_by(altpath, b) %>%
		summarise(p = sum(p))
	res <- rnorm(nboot, gx, gx.se) - altpath$p
	pe <- gx - sum(gp * px)
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


