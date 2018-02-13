

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

