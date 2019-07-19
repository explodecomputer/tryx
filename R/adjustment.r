

#' Outlier adjustment estimation
#' 
#' How much of the heterogeneity due to the outlier can be explained by alternative pathways?
#' 
#' @param tryxscan Output from \code{tryx.scan}
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#'
#' @export
#' @return data frame of adjusted effect estimates and heterogeneity stats
tryx.adjustment <- function(tryxscan, id_remove=NULL)
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
	sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
	sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)
	sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)

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
				message("x<-p->y:  ", a$SNP, "\t- ", sig$outcome[i])
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
			}
			temp <- dat
			temp$beta.exposure[temp$SNP == a$SNP] <- a$adj.beta.exposure
			temp$se.exposure[temp$SNP == a$SNP] <- a$adj.se.exposure
			temp$beta.outcome[temp$SNP == a$SNP] <- a$adj.beta.outcome
			temp$se.outcome[temp$SNP == a$SNP] <- a$adj.se.outcome
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


tryx.adjustment.mv <- function(tryxscan, id_remove=NULL, proxies=FALSE)
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
	sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
	sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)
	sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)

	id.exposure <- tryxscan$dat$id.exposure[1]
	id.outcome <- tryxscan$dat$id.outcome[1]

#	for each outlier SNP find its effects on 

	sigo1 <- subset(sig, id.outcome %in% sigo$id.exposure) %>% group_by(SNP) %>% mutate(snpcount=n()) %>% arrange(desc(snpcount), SNP)
	snplist <- unique(sigo1$SNP)
	mvo <- list()
	dat <- tryxscan$dat
	dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
	dat$Q <- sum(dat$qi)
	dat$orig.beta.outcome <- dat$beta.outcome
	dat$orig.se.outcome <- dat$se.outcome
	dat$orig.qi <- dat$qi
	dat$outlier <- FALSE
	dat$outlier[dat$SNP %in% tryxscan$outliers] <- TRUE
	for(i in 1:length(snplist))
	{
		message("Estimating joint effects of the following traits associated with ", snplist[i])
		temp <- subset(sigo1, SNP %in% snplist[i])
		message(paste(unique(temp$outcome), collapse="\n"))
		mvexp <- suppressMessages(mv_extract_exposures(c(id.exposure, temp$id.outcome), find_proxies=proxies))
		mvout <- suppressMessages(extract_outcome_data(mvexp$SNP, id.outcome))
		mvdat <- suppressMessages(mv_harmonise_data(mvexp, mvout))
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
		mvo[[i]] <- mv_multiple(mvdat2)$result
		mvo[[i]] <- subset(mvo[[i]], !exposure %in% tryxscan$dat$exposure[1])
		temp2 <- with(temp, tibble(SNP=SNP, exposure=outcome, snpeff=beta.outcome, snpeff.se=se.outcome, snpeff.pval=pval.outcome))
		mvo[[i]] <- merge(mvo[[i]], temp2, by="exposure")
		boo <- with(subset(dat, SNP == snplist[i]), 
			bootstrap_path(
				beta.outcome,
				se.outcome,
				mvo[[i]]$snpeff,
				mvo[[i]]$snpeff.se,
				mvo[[i]]$b,
				mvo[[i]]$se
			))
		dat$beta.outcome[dat$SNP == snplist[i]] <- boo[1]
		dat$se.outcome[dat$SNP == snplist[i]] <- boo[2]
	}
	mvo <- bind_rows(mvo)

	dat$qi[dat$mr_keep] <- cochrans_q(dat$beta.outcome[dat$mr_keep] / dat$beta.exposure[dat$mr_keep], dat$se.outcome[dat$mr_keep] / abs(dat$beta.exposure[dat$mr_keep]))
	return(list(mvo=mvo, dat=dat))
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
#' @param duplicate_outliers_method Sometimes more than one trait will associate with a particular outlier. TRUE = only keep the trait that has the biggest influence on heterogeneity
#' 
#' @export
#' @return List of 
#' - adj_full: data frame of SNP adjustments for all candidate traits
#' - adj: The results from adj_full selected to adjust the exposure-outcome model
#' - Q: Heterogeneity stats
#' - estimates: Adjusted and unadjested exposure-outcome effects
#' - plot: Radial plot showing the comparison of different methods and the changes in SNP effects ater adjustment
tryx.analyse <- function(tryxscan, plot=TRUE, id_remove=NULL, filter_duplicate_outliers=TRUE)
{

	analysis <- list()
	adj_full <- tryx.adjustment(tryxscan, id_remove)
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
	# if("simulation" %in% names(tryxscan))
	# {
	# 	detection <- list()
	# 	detection$nu1_correct <- sum(adj$SNP %in% 1:tryxscan$simulation$nu1 & adj$what != "p->y")
	# 	detection$nu1_incorrect <- sum(adj$SNP %in% 1:tryxscan$simulation$nu1 & adj$what == "p->y")
	# 	detection$nu2_correct <- sum(adj$SNP %in% (1:tryxscan$simulation$nu2 + tryxscan$simulation$nu1) & adj$what == "p->y")
	# 	detection$nu2_incorrect <- sum(adj$SNP %in% (1:tryxscan$simulation$nu2 + tryxscan$simulation$nu1) & adj$what != "p->y")
	# 	detection$no_outlier_flag <- tryxscan$simulation$no_outlier_flag
	# 	analysis$detection <- detection
	# }

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

	dat_rem <- subset(dat, !SNP %in% temp$SNP)
	dat_rem2 <- subset(dat, !SNP %in% tryxscan$outliers)
	# adjust everything that's possible remove remaining outliers
	dat_rem3 <- rbind(temp, dat_rem2) %>% filter(!duplicated(SNP))



	estimates <- tibble()

	
	# Raw IVW
	tt <- dat
	mod <- summary(lm(ratiow ~ -1 + weights, data=tt))
	estimates <- bind_rows(estimates, 
		tibble(
			est="Raw",
			b=coefficients(mod)[1,1], 
			se = coefficients(mod)[1,2],
			pval=coefficients(mod)[1,4],
			nsnp = nrow(tt),
			Q = sum(tt$qi),
			int = 0
		)
	)

	# Outliers removed (all)
	tt <- subset(dat, !SNP %in% tryxscan$outliers)
	mod <- try(summary(lm(ratiow ~ -1 + weights, data=tt)))
	if(class(mod) != "try-error")
	{
		estimates <- bind_rows(estimates, 
			tibble(
				est="Outliers removed (all)",
				b=coefficients(mod)[1,1], 
				se = coefficients(mod)[1,2],
				pval=coefficients(mod)[1,4],
				nsnp = nrow(tt),
				Q = sum(tt$qi),
				int = 0
			)
		)
	}

	# Outliers removed (candidates)
	tt <- subset(dat, !SNP %in% temp$SNP)
	mod <- try(summary(lm(ratiow ~ -1 + weights, data=tt)))
	if(class(mod) != "try-error")
	{
		estimates <- bind_rows(estimates, 
			tibble(
				est="Outliers removed (candidates)",
				b=coefficients(mod)[1,1], 
				se = coefficients(mod)[1,2],
				pval=coefficients(mod)[1,4],
				nsnp = nrow(tt),
				Q = sum(tt$qi),
				int = 0
			)
		)
	}

	# Outliers adjusted
	tt <- rbind(temp, dat) %>% filter(!duplicated(SNP))
	tt$qi <- cochrans_q(tt$beta.outcome / tt$beta.exposure, tt$se.outcome / abs(tt$beta.exposure))
	analysis$Q$adj_Q <- sum(tt$qi)
	mod <- try(summary(lm(ratiow ~ -1 + weights, data=tt)))
	if(class(mod) != "try-error")
	{
		estimates <- bind_rows(estimates, 
			tibble(
				est="Outliers adjusted",
				b=coefficients(mod)[1,1], 
				se = coefficients(mod)[1,2],
				pval=coefficients(mod)[1,4],
				nsnp = nrow(tt),
				Q = sum(tt$qi),
				int = 0
			)
		)
	}

	if("mvres" %in% names(tryxscan))
	{
		estimates <- bind_rows(estimates, 
			tibble(
				est="Multivariable MR",
				b=tryxscan$mvres$result$b[1], 
				se = tryxscan$mvres$result$se[1],
				pval = tryxscan$mvres$result$pval[1],
				nsnp = tryxscan$mvres$result$nsnp[1],
				Q = NA,
				int = 0
			)
		)
	}

	if("true_outliers" %in% names(tryxscan))
	{
		tt <- subset(dat, !SNP %in% tryxscan$true_outliers)
		mod <- try(summary(lm(ratiow ~ -1 + weights, data=tt)))
		if(class(mod) != "try-error")
		{
			estimates <- bind_rows(estimates,
				tibble(
					est=c("Oracle"),
					b=c(coefficients(mod)[1,1]), 
					se=c(coefficients(mod)[1,2]), 
					pval=c(coefficients(mod)[1,4]),
					nsnp=c(nrow(tt)),
					Q = c(sum(tt$qi)),
					int=0
				)
			)
		}
	}


	estimates <- bind_rows(estimates)
	estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q) 

	analysis$estimates <- estimates


	# est3 <- summary(lm(ratiow ~ -1 + weights, data=dat_adj))
	# est4 <- try(summary(lm(ratiow ~ -1 + weights, data=dat_rem2)))
	# est5 <- try(summary(lm(ratiow ~ -1 + weights, data=dat_rem3)))
	# if(class(est4) == "try-error")
	# {
	# 	est4 <- list(coefficients = matrix(NA, 2,4))
	# }

	# estimates <- data.frame(
	# 	est=c("Raw", "Outliers removed (candidates)", "Outliers removed (all)", "Outliers adjusted", "Mixed", "Multivariable MR"),
	# 	b=c(coefficients(est2)[1,1], coefficients(est1)[1,1], coefficients(est4)[1,1], coefficients(est3)[1,1], coefficients(est5)[1,1], tryxscan$mvres$result$b[1]), 
	# 	se=c(coefficients(est2)[1,2], coefficients(est1)[1,2], coefficients(est4)[1,2], coefficients(est3)[1,2], coefficients(est5)[1,2], tryxscan$mvres$result$se[1]), 
	# 	pval=c(coefficients(est2)[1,4], coefficients(est1)[1,4], coefficients(est4)[1,4], coefficients(est3)[1,4], coefficients(est5)[1,4], tryxscan$mvres$result$pval[1]),
	# 	nsnp=c(nrow(dat), nrow(dat_rem), nrow(dat_rem2), nrow(dat_adj), nrow(dat_rem3), tryxscan$mvres$result$nsnp[1]),
	# 	Q = c(sum(dat$qi), sum(dat_rem$qi), sum(dat_rem2$qi), sum(dat_adj$qi), sum(dat_rem3$qi), NA),
	# 	int=0
	# )
	# estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q) 

	# estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q) 

	# analysis$estimates <- estimates

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



#' Adjust and analyse the tryx results
#'
#' Similar to tryx.analyse, but when there are multiple traits associated with a single variant then we use a LASSO-based multivariable approach 
#'
#' @param tryxscan Output from \code{tryx.scan}
#' @param plot Whether to plot or not. Default is TRUE
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#' @param proxies Look for proxies in the MVMR methods. Default = FALSE.
#'
#' @export
#' @return List of 
#' - adj_full: data frame of SNP adjustments for all candidate traits
#' - adj: The results from adj_full selected to adjust the exposure-outcome model
#' - Q: Heterogeneity stats
#' - estimates: Adjusted and unadjested exposure-outcome effects
#' - plot: Radial plot showing the comparison of different methods and the changes in SNP effects ater adjustment
tryx.analyse.mv <- function(tryxscan, plot=TRUE, id_remove=NULL, proxies=FALSE)
{
	adj <- tryx.adjustment.mv(tryxscan, id_remove, proxies)
	dat <- subset(adj$dat, mr_keep)
	dat$orig.ratio <- dat$orig.beta.outcome / dat$beta.exposure
	dat$orig.weights <- sqrt(dat$beta.exposure^2 / dat$orig.se.outcome^2)
	dat$orig.ratiow <- dat$orig.ratio * dat$orig.weights
	dat$ratio <- dat$beta.outcome / dat$beta.exposure
	dat$weights <- sqrt(dat$beta.exposure^2 / dat$se.outcome^2)
	dat$ratiow <- dat$ratio * dat$weights

	est_raw <- summary(lm(orig.ratiow ~ -1 + orig.weights, data=dat))
	est_adj <- summary(lm(ratiow ~ -1 + weights, data=dat))
	est_out1 <- summary(lm(orig.ratiow ~ -1 + orig.weights, data=subset(dat, !SNP %in% tryxscan$outliers)))
	est_out2 <- summary(lm(orig.ratiow ~ -1 + orig.weights, data=subset(dat, !SNP %in% adj$mvo$SNP)))

	estimates <- data.frame(
		est=c("Raw", "Outliers removed (all)", "Outliers removed (candidates)", "Outliers adjusted"),
		b=c(coefficients(est_raw)[1,1], coefficients(est_out2)[1,1], coefficients(est_out1)[1,1], coefficients(est_adj)[1,1]), 
		se=c(coefficients(est_raw)[1,2], coefficients(est_out2)[1,2], coefficients(est_out1)[1,2], coefficients(est_adj)[1,2]), 
		pval=c(coefficients(est_raw)[1,4], coefficients(est_out2)[1,4], coefficients(est_out1)[1,4], coefficients(est_adj)[1,4]),
		nsnp=c(nrow(dat), nrow(subset(dat, !SNP %in% tryxscan$outliers)), nrow(subset(dat, !SNP %in% adj$mvo$SNP)), nrow(dat)),
		Q = c(sum(dat$orig.qi), sum(subset(dat, !SNP %in% tryxscan$outliers)$qi), sum(subset(dat, !SNP %in% adj$mvo$SNP)$qi), sum(dat$qi)),
		int=0
	)
	estimates$Isq <- pmax(0, (estimates$Q - estimates$nsnp - 1) / estimates$Q) 

	analysis <- list(
		estimates=estimates,
		mvo=adj$mvo,
		dat=adj$dat
	)

	if(plot)
	{	
		dato <- subset(dat, SNP %in% tryxscan$outliers)
		datadj <- subset(dat, SNP %in% adj$mvo$SNP)
		datadj$x <- datadj$orig.weights
		datadj$xend <- datadj$weights
		datadj$y <- datadj$orig.ratiow
		datadj$yend <- datadj$ratiow

		mvog <- tidyr::separate(adj$mvo, exposure, sep="\\|", c("exposure", "temp", "id"))
		mvog$exposure <- gsub(" $", "", mvog$exposure)
		mvog <- group_by(mvog, SNP) %>% summarise(label = paste(exposure, collapse="\n"))
		datadj <- merge(datadj, mvog, by="SNP")

		labs <- rbind(
			tibble(label=dato$SNP, x=dato$orig.weights, y=dato$orig.ratiow, col="grey50"),
			tibble(label=datadj$label, x=datadj$weights, y=datadj$ratiow, col="grey100")
		)


		p <- ggplot(dat, aes(y=orig.ratiow, x=orig.weights)) +
		geom_abline(data=estimates, aes(slope=b, intercept=0, linetype=est)) +
		geom_label_repel(data=labs, aes(label=label, x=x, y=y, colour=col), size=2, segment.color = "grey50", show.legend = FALSE) +
		geom_point(data=dato, size=4) +
		# geom_label_repel(data=labs, aes(label=label, x=weights, y=ratiow), size=2, segment.color = "grey50") +
		geom_point(data=datadj, aes(x=weights, y=ratiow)) +
		geom_point(aes(colour=SNP %in% dato$SNP), show.legend=FALSE) +
		geom_segment(data=datadj, colour="grey50", aes(x=x, xend=xend, y=y, yend=yend), arrow = arrow(length = unit(0.01, "npc"))) +
		labs(colour=NULL, linetype="Estimate", x="w", y="beta * w")
		# xlim(c(0, max(dat$weights))) +
		# ylim(c(min(0, dat$ratiow, temp$ratiow), max(dat$ratiow, temp$ratiow)))
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
	altpath <- tibble(
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
