#' Simulate data to test tryx
#' 
#' Choose a sample size, number of horizontal pleiotropy paths, and number of confounders. Create summary data for all necessary traits, then perform tryx scan
#' 
#' @param nid Sample size to simulate
#' @param nu1 Number of traits which influence both X and Y
#' @param nu2 Number of traits which influence Y
#' @param bxy Causal effect of X on Y
#' @param outliers_known Assume that all outliers are detected (defaukt = TRUE). If FALSE then Radial MR is used to detect outliers
#' 
#' @export
#' @return Results from simulations and tryx scan
simulate.tryx <- function(nid, nu1, nu2, bxy=3, outliers_known=TRUE)
{
	# scenario 1 - confounder g -> u; u -> x; u -> y
	# scenario 2 - pleiotropy g -> u -> y

	out <- list()

	cpg <- require(simulateGP)
	if(!cpg)
	{
		stop("Please install the simulateGP package\ndevtools::install_github('explodecomputer/simulateGP')")
	}

	ngx <- 30
	ngu1 <- 30
	ngu2 <- 30

	gx <- make_geno(nid, ngx, 0.5)
	nu <- nu1 + nu2

	# Effect sizes
	bgu1 <- lapply(1:nu1, function(x) runif(ngu1))
	bgu2 <- lapply(1:nu2, function(x) runif(ngu2))
	bgx <- runif(ngx)
	bu1x <- rnorm(nu1)
	bu1y <- rnorm(nu1)
	bu2y <- rnorm(nu2)
	bxy <- bxy

	gu1 <- list()
	u1 <- matrix(nid * nu1, nid, nu1)
	for(i in 1:nu1)
	{
		gu1[[i]] <- cbind(gx[,i], make_geno(nid, ngu1-1, 0.5))
		u1[,i] <- gu1[[i]] %*% bgu1[[i]] + rnorm(nid)
	}
	gu2 <- list()
	u2 <- matrix(nid * nu2, nid, nu2)
	for(j in 1:nu2)
	{
		gu2[[j]] <- cbind(gx[,i+j], make_geno(nid, ngu1-1, 0.5))
		u2[,j] <- gu2[[j]] %*% bgu2[[j]] + rnorm(nid)
	}

	ux <- rep(2, nu1)
	x <- drop(gx %*% bgx + u1 %*% bu1x + rnorm(nid))
	# u3 <- gu %*% c(0.5, runif(6)) + x * -3 + rnorm(nid, sd=6)
	y <- drop(x * bxy + u1 %*% bu1y + u2 %*% bu2y + rnorm(nid))

	out$dat <- get_effs(x, y, gx, "X", "Y")
	out$dat$mr_keep <- TRUE
	no_outlier_flag <- FALSE

	if(outliers_known)
	{
		outliers <- as.character(c(1:nu))
	} else {
		radial <- RadialMR::ivw_radial(RadialMR::format_radial(out$dat$beta.exposure, out$dat$beta.outcome, out$dat$se.exposure, out$dat$se.outcome, out$dat$SNP), 0.05/nrow(out$dat), summary=FALSE)
		if(radial$outliers[1] == "No significant outliers")
		{
			outliers <- as.character(c(1:nu))
			no_outlier_flag <- TRUE
		} else {
			outliers <- as.character(sort(radial$outliers$SNP))
			message("Number of outliers: ", length(outliers))
		}
	}
	# output$radialmr <- radial
	nout <- length(outliers)
	outliers <- subset(out$dat, SNP %in% outliers)$SNP
	out$outliers <- outliers
	out$id_list <- c(paste0("U1_", 1:nu1), paste0("U2_", 1:nu2))

	# temp1 <- gwas(u1, gx[,1:2])

	temp <- list()
	for(i in 1:nu1)
	{
		temp[[i]] <- gwas(u1[,i], gx[,outliers, drop=FALSE])
		temp[[i]]$SNP <- outliers
		temp[[i]]$outcome <- paste0("U1_", i)
	}
	for(i in 1:nu2)
	{
		temp[[i + nu1]] <- gwas(u2[,i], gx[,outliers, drop=FALSE])
		temp[[i + nu1]]$SNP <- outliers
		temp[[i + nu1]]$outcome <- paste0("U2_", i)
	}

	temp <- bind_rows(temp)

	out$search <- data.frame(SNP=temp$SNP, beta.outcome=temp$bhat, se.outcome=temp$se, pval.outcome=temp$pval, id.outcome=temp$outcome, outcome=temp$outcome, mr_keep=TRUE, effect_allele.outcome="A", other_allele.outcome="G", eaf.outcome=0.5)
	out$search$sig <- out$search$pval.outcome < 5e-8


	## get effects
	u1x <- list()
	u1y <- list()
	for(i in 1:nu1)
	{
		u1x[[i]] <- get_effs(u1[,i], x, gu1[[i]], paste0("U1_", i), "X")
		u1y[[i]] <- get_effs(u1[,i], y, gu1[[i]], paste0("U1_", i), "Y")
		u1x[[i]]$SNP <- u1y[[i]]$SNP <- c(i, paste0("u1_",1:(ngu1-1)))
		u1x[[i]]$effect_allele.exposure <- u1y[[i]]$effect_allele.exposure <- "A"
		u1x[[i]]$other_allele.exposure <- u1y[[i]]$other_allele.exposure <- "G"
		u1x[[i]]$effect_allele.outcome <- u1y[[i]]$effect_allele.outcome <- "A"
		u1x[[i]]$other_allele.outcome <- u1y[[i]]$other_allele.outcome <- "G"
		u1x[[i]]$eaf.exposure <- u1y[[i]]$eaf.exposure <- 0.5
		u1x[[i]]$eaf.outcome <- u1y[[i]]$eaf.outcome <- 0.5
	}

	u2x <- list()
	u2y <- list()
	for(i in 1:nu2)
	{
		u2x[[i]] <- get_effs(u2[,i], x, gu2[[i]], paste0("U2_", i), "X")
		u2y[[i]] <- get_effs(u2[,i], y, gu2[[i]], paste0("U2_", i), "Y")
		u2x[[i]]$SNP <- u2y[[i]]$SNP <- c(i+nu1, paste0("u2_",1:(ngu2-1)))
		u2x[[i]]$effect_allele.exposure <- u2y[[i]]$effect_allele.exposure <- "A"
		u2x[[i]]$other_allele.exposure <- u2y[[i]]$other_allele.exposure <- "G"
		u2x[[i]]$effect_allele.outcome <- u2y[[i]]$effect_allele.outcome <- "A"
		u2x[[i]]$other_allele.outcome <- u2y[[i]]$other_allele.outcome <- "G"
		u2x[[i]]$eaf.exposure <- u2y[[i]]$eaf.exposure <- 0.5
		u2x[[i]]$eaf.outcome <- u2y[[i]]$eaf.outcome <- 0.5
	}



	out$candidate_instruments <- subset(rbind(bind_rows(u1x), bind_rows(u2x)), select=c(
		SNP, beta.exposure, se.exposure, id.exposure, exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, pval.exposure
	))
	out2 <- subset(out$search, sig)
	out$candidate_instruments <- group_by(out$candidate_instruments, id.exposure) %>%
		do({
			x <- .
			y <- subset(out2, id.outcome == x$id.exposure[1])
			x <- subset(x, !SNP %in% y$SNP)
			x
		})

	out$candidate_outcome <- subset(rbind(bind_rows(u1y), bind_rows(u2y)), select=c(
		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
	))

	out$candidate_outcome_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_outcome))
	out$candidate_outcome_mr <- suppressMessages(mr(out$candidate_outcome_dat, method_list="mr_ivw"))


	out$candidate_exposure <- subset(rbind(bind_rows(u1x), bind_rows(u2x)), select=c(
		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
	))

	out$candidate_exposure_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_exposure))
	out$candidate_exposure_mr <- suppressMessages(mr(out$candidate_exposure_dat, method_list="mr_ivw"))

	out$simulation <- list()
	out$simulation$no_outlier_flag <- no_outlier_flag
	out$simulation$nu1 <- nu1
	out$simulation$nu2 <- nu2
	out$simulation$outliers_known <- outliers_known
	out$simulation$bxy <- bxy

	return(out)
}

