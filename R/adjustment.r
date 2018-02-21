

#' Outlier adjustment estimation
#' 
#' How much of the heterogeneity due to the outlier can be explained by alternative pathways?
#' 
#' @param tryxscan Output from \code{tryx.scan}
#' @export
#' @return data frame of adjusted effect estimates and heterogeneity stats
tryx.adjustment <- function(tryxscan)
{
	# for each outlier find the candidate MR analyses
	# if only exposure then ignore
	# if only outcome then re-estimate the snp-outcome association
	# if exposure and outcome then re-estimate the snp-exposure and snp-outcome association
	# outlier
	# candidate
	# what 
	# old.beta.exposure
	# adj.beta.exposure
	# old.beta.outcome
	# adj.beta.outcome
	# candidate.beta.outcome
	# candidate.se.outcome
	# candidate.beta.exposure
	# candidate.se.exposure
	# old.deviation
	# new.deviation

	if(!any(tryxscan$search$sig))
	{
		return(tibble())
	}

	l <- list()
	sig <- subset(tryxscan$search, sig)
	sige <- subset(tryxscan$candidate_exposure_mr, sig)
	sigo <- subset(tryxscan$candidate_outcome_mr, sig)

	dat <- tryxscan$dat
	dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
	dat$Q <- sum(dat$qi)

	for(i in 1:nrow(sig))
	{
		a <- subset(dat, SNP == sig$SNP[i], select=c(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome, qi, Q))
		if(sig$id.outcome[i] %in% sigo$id.exposure)
		{
			a$candidate <- sig$outcome[i]
			a$i <- i
			if(sig$id.outcome[i] %in% sige$id.exposure) {
				message("p->x; p->y: ", a$SNP, " - ", sig$outcome[i])
				a$what <- "p->x; p->y"
				a$candidate.beta.exposure <- sige$b[sige$id.exposure == sig$id.outcome[i]]
				a$candidate.se.exposure <- sige$se[sige$id.exposure == sig$id.outcome[i]]
				a$candidate.beta.outcome <- sigo$b[sigo$id.exposure == sig$id.outcome[i]]
				a$candidate.se.outcome <- sigo$se[sigo$id.exposure == sig$id.outcome[i]]
				b <- bootstrap_path(a$beta.exposure, a$se.exposure, sig$beta.outcome[i], sig$se.outcome[i], sige$b[sige$id.exposure == sig$id.outcome[i]], sige$se[sige$id.exposure == sig$id.outcome[i]])
				a$adj.beta.exposure <- b[1]
				a$adj.se.exposure <- b[2]
				b <- bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
				a$adj.beta.outcome <- b[1]
				a$adj.se.outcome <- b[2]
			} else {
				message("p->y: ", a$SNP, " - ", sig$outcome[i])
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
			}
			temp <- dat
			temp$beta.exposure[temp$SNP == a$SNP] <- a$adj.beta.exposure
			temp$beta.exposure.se[temp$SNP == a$SNP] <- a$adj.beta.exposure.se
			temp$beta.outcome[temp$SNP == a$SNP] <- a$adj.beta.outcome
			temp$beta.outcome.se[temp$SNP == a$SNP] <- a$adj.beta.outcome.se
			temp$qi <- cochrans_q(temp$beta.outcome / temp$beta.exposure, temp$se.outcome / abs(temp$beta.exposure))
			a$adj.qi <- temp$qi[temp$SNP == a$SNP]
			a$adj.Q <- sum(temp$qi)

			l[[i]] <- a
		}

	}
	l <- bind_rows(l)
	l$d <- (l$adj.qi / l$adj.Q) / (l$qi / l$Q)
	return(l)
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
#' 
#' @export
#' @return List of 
#' - adj_full: data frame of SNP adjustments for all candidate traits
#' - adj: The results from adj_full selected to adjust the exposure-outcome model
#' - Q: Heterogeneity stats
#' - estimates: Adjusted and unadjested exposure-outcome effects
#' - plot: Radial plot showing the comparison of different methods and the changes in SNP effects ater adjustment
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
	dat_rem2 <- subset(dat, !SNP %in% tryxscan$outliers)

	est1 <- summary(lm(ratiow ~ -1 + weights, data=dat_rem))
	est2 <- summary(lm(ratiow ~ -1 + weights, data=dat))
	est3 <- summary(lm(ratiow ~ -1 + weights, data=dat_adj))
	est4 <- summary(lm(ratiow ~ -1 + weights, data=dat_rem2))

	estimates <- data.frame(
		est=c("Raw", "Outliers removed (all)", "Outliers removed (candidates)", "Outliers adjusted"),
		b=c(coefficients(est2)[1,1], coefficients(est1)[1,1], coefficients(est4)[1,1], coefficients(est3)[1,1]), 
		se=c(coefficients(est2)[1,2], coefficients(est1)[1,2], coefficients(est4)[1,2], coefficients(est3)[1,2]), 
		pval=c(coefficients(est2)[1,4], coefficients(est1)[1,4], coefficients(est4)[1,4], coefficients(est3)[1,4]),
		nsnp=c(nrow(dat), nrow(dat_rem), nrow(dat_rem2), nrow(dat_adj)),
		Q = c(sum(dat$qi), sum(dat_rem$qi), sum(dat_rem2$qi), sum(dat_adj$qi)),
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
	res <- rnorm(nboot, gx, gx.se) - rnorm(nboot, gp, gp.se) * rnorm(nboot, px, px.se)
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
