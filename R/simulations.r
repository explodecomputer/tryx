#' Simulate data to test tryx
#' 
#' # Create data using model here: https://photos.google.com/photo/AF1QipNbxPRZyr7iODKRwuwrEmYALTpxnyudq7xovC7v. Create summary data for all necessary traits, then perform tryx scan
#' x = exposure
#' y = outcome
#' u1 = confounder of x and y
#' u2 = mediator of horizontal pleiotropy from gx to y. One for each gx
#' u3 = mediator between x and y
#'
#' @param nid = 10000 Number of samples
#' @param ngx = 30 Number of direct instruments to x
#' @param ngu1 = 30 Number of instruments influencing confounder of x and y
#' @param ngu2 = 30 Number of instruments per mediating
#' @param nu2 = 2 Number of gx that are pleiotropic
#' @param ngu3 = 30 Number of instruments for u3 mediator
#' @param vgx = 0.2 Variance explained by gx instruments
#' @param vgu1 = 0.6 Variance explained by u1 instruments
#' @param vgu2 = 0.2 Variance explained by u2 instruments
#' @param vgu3 = 0.2 Variance explained by u3 instruments
#' @param bxy = 0 Causal effect of x on y
#' @param bu1x = 0.6 Effect of u1 on x
#' @param bu1y = 0.4 Effect of u1 on y
#' @param bxu3 = 0.3 Effect of x on u3
#' @param bu3y = 0 Effect of u3 on y
#' @param vgxu2 = 0.2 Variance explained by each gx instrument on each u2 mediator
#' @param vu2y = 0.2 Variance explained by all u2 mediators on y
#' @param mininum_instruments = 10 Minimum number of instruments required to have been detected to run simulation
#' @param instrument_threshold = "bonferroni" Threshold, either numeric or 'bonferroni'
#' @param outlier_threshold = "bonferroni" Threshold, either numeric or 'bonferroni'
#' @param outliers_known = "detected" Either detected = using radial mr to find heterogeneity outliers; known = adjust for all known invalid instruments; all = adjust for all instruments
#' @param directional_bias = FALSE Is the pleiotropic effect randomly centered around 0 (FALSE) or does it have a non-0 mean (TRUE)
#'
#' @export
#' @return list for tryx.analyse
tryx.simulate <- function(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 2, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0, bu1x = 0.6, bu1y = 0.4, bxu3 = 0.3, bu3y = 0, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = "bonferroni", outliers_known = "detected", directional_bias = FALSE)
{
	out <- list()

	message("Generating genetic effects")

	gx <- make_geno(nid, ngx, 0.5)
	gu1 <- make_geno(nid, ngu1, 0.5)
	gu2 <- make_geno(nid, ngu2 * nu2, 0.5)
	gu3 <- make_geno(nid, ngu3, 0.5)
	u1 <- make_phen(choose_effects(ngu1, vgu1), gu1)
	x <- make_phen(c(abs(choose_effects(ngx, vgx)), bu1x), cbind(gx, u1))
	u3 <- make_phen(c(choose_effects(ngu3, vgu3), bxu3), cbind(gu3, x))


	# Some number of the direct gx SNPs also influence the outcome via a mediator
	# This is mediated pleiotropy

	message("Creating potential pleiotropy paths")

	u2 <- matrix(rnorm(nid * nu2), nid, nu2)
	for(i in 1:nu2)
	{
		message("SNP ", i, " is pleiotropic")
		u2[,i] <- make_phen(
			c(choose_effects(ngu2, vgu2), abs(choose_effects(1, vgxu2))),
			cbind(gu2[,((i-1) * ngu2 + 1):(ngu2 * i)], gx[,i])
		)
	}

	message("Creating outcome")

	bu2y <- choose_effects(nu2, vu2y)

	if(directional_bias) 
	{
		message("U2 bias is directional")
		bu2y <- abs(bu2y)
	}
	y <- make_phen(
		c(bxy, bu3y, bu1y, bu2y),
		cbind(x, u3, u1, u2)
	)

	# Single genotype data source for simplicity

	G <- cbind(gx, gu1, gu2, gu3)
	G_annotation <- tibble(
		snp=1:ncol(G),
		origin1=c(
			rep("x", ngx),
			rep("u1", ngu1),
			rep(paste0("u2_",1:nu2), each=ngu2),
			rep("u3", ngu3)
		),
		origin2=c(
			rep("x", ngx),
			rep("u1", ngu1),
			rep("u2", each=ngu2*nu2),
			rep("u3", ngu3)
		)
	)

	# How do we select instruments?

	if(instrument_threshold == "bonferroni")
	{
		instrument_threshold <- 0.05 / ncol(G)
	} else {
		stopifnot(is.numeric(instrument_threshold))
	}
	message("instrument_threshold = ", instrument_threshold)


	# Get summary set

	message("Getting instruments")

	out$dat_all <- get_effs(x, y, G, "X", "Y")
	out$dat <- subset(out$dat_all, pval.exposure < instrument_threshold)
	stopifnot(nrow(out$dat) > mininum_instruments)
	out$dat$mr_keep <- TRUE
	G_annotation$detected <- G_annotation$snp %in% out$dat$SNP
	message("Detected ", nrow(out$dat), " associations for x")

	## Determine valid instruments

	if(vu2y == 0)
	{
		valid <- 1:ngx
	} else {
		valid <- (1:ngx)[! (1:ngx %in% 1:nu2)]
	}

	if(bu1y == 0)
	{
		valid <- c(valid, (ngx+1):(ngx+ngu1))
	}
	message("Valid instruments: ", paste(valid, collapse=","))
	G_annotation$valid <- G_annotation$snp %in% valid

	## known invalid instruments
	invalid <- which(! 1:ncol(G) %in% valid)

	out$valid <- valid
	out$invalid <- invalid



	# Get outliers

	no_outlier_flag <- FALSE
	if(outliers_known == "known")
	{
		message("Assuming knowledge of all outliers")
		outliers <- out$dat$SNP[out$dat$SNP %in% invalid]
	} else if(outliers_known %in% c("detected", "all")) {
		message("Detecting outliers")
		radial <- RadialMR::ivw_radial(RadialMR::format_radial(out$dat$beta.exposure, out$dat$beta.outcome, out$dat$se.exposure, out$dat$se.outcome, out$dat$SNP), ifelse(outlier_threshold == "bonferroni", 0.05/nrow(out$dat), outlier_threshold), weights=3)
		if(radial$outliers[1] == "No significant outliers")
		{
			message("No significant outliers detected")
			outliers <- as.character(out$dat$SNP[out$dat$SNP %in% invalid])
			if(length(outliers) == 0) outliers <- 1
			no_outlier_flag <- TRUE
		} else {
			outliers <- as.character(sort(radial$outliers$SNP))
			message("Number of outliers: ", length(outliers))
		}
	} else {
		stop("outliers_known must be 'known', 'detected' or 'all'")
	}

	nout <- length(outliers)
	outliers <- subset(out$dat, SNP %in% outliers)$SNP
	out$outliers <- outliers
	if(outliers_known == "all")
	{
		message("Treating all variants as outliers")
		outliers <- out$dat$SNP
	}
	# out$id_list <- c(paste0("U1_", 1:nu1), paste0("U2_", 1:nu2))
	message("Outliers to be analysed: ", paste(outliers, collapse=","))

	message("Perform outlier scan")

	U <- cbind(u1, u2, u3)
	colnames(U) <- c("U1", paste0("U2_", 1:ncol(u2)), "U3")
	outlier_scan <- list()
	for(i in 1:ncol(U))
	{
		outlier_scan[[colnames(U)[i]]] <- gwas(U[,i], G[,outliers, drop=FALSE])
		outlier_scan[[colnames(U)[i]]]$SNP <- outliers
		outlier_scan[[colnames(U)[i]]]$outcome <- colnames(U)[i]
	}
	outlier_scan <- bind_rows(outlier_scan)
	out$search <- tibble::tibble(
		SNP=outlier_scan$SNP, 
		beta.outcome=outlier_scan$bhat, 
		se.outcome=outlier_scan$se, 
		pval.outcome=outlier_scan$pval, 
		id.outcome=outlier_scan$outcome, 
		outcome=outlier_scan$outcome, 
		mr_keep=TRUE, 
		effect_allele.outcome="A", 
		other_allele.outcome="G", 
		eaf.outcome=0.5
	)
	out$search$sig <- out$search$pval.outcome < 5e-8
	message("Found the following candidate traits: ", paste(unique(subset(out$search, sig)$id.outcome), collapse=","))


	message("Organise data for mvmr")

	traitlist <- as.character(unique(subset(out$search, sig)$id.outcome))
	PHEN <- cbind(x, U[,traitlist])
	PHEN <- lapply(seq_len(ncol(PHEN)), function(i) PHEN[,i])
	names(PHEN) <- c("x", traitlist)
	message("Analysing ", length(PHEN), " traits: ", paste(names(PHEN), collapse=","))

	snplist <- sapply(1:length(PHEN), function(x) {
		subset(gwas(PHEN[[x]], G), pval < instrument_threshold)$snp
	}) %>% unlist() %>% sort %>% unique
	message("Found ", length(snplist), " unique instruments")

	message("Perform mvmr")
	out$mvres <- try({
		mvdat <- make_mvdat(PHEN, y, G[,snplist])
		mvres <- mv_multiple(mvdat)
		mvres$result$exposure <- c("x", traitlist)
		mvres
	})


	message("Building effects for ", length(traitlist), " trait(s)")

	## get effects
	ux <- list()
	uy <- list()
	for(i in traitlist)
	{
		message(i)
		ux[[i]] <- get_effs(U[,i], x, G, i, "X") %>% subset(pval.exposure < instrument_threshold)
		uy[[i]] <- get_effs(U[,i], y, G, i, "Y") %>% subset(pval.exposure < instrument_threshold)
		if(nrow(ux[[i]]) > 0)
		{
			ux[[i]]$effect_allele.exposure <- "A"
			ux[[i]]$other_allele.exposure <- "G"
			ux[[i]]$effect_allele.outcome <- "A"
			ux[[i]]$other_allele.outcome <- "G"
			ux[[i]]$eaf.exposure <- 0.5
			ux[[i]]$eaf.outcome <- 0.5
		}

		if(nrow(uy[[i]]) > 0)
		{
			uy[[i]]$effect_allele.exposure <- "A"
			uy[[i]]$other_allele.exposure <- "G"
			uy[[i]]$effect_allele.outcome <- "A"
			uy[[i]]$other_allele.outcome <- "G"
			uy[[i]]$eaf.exposure <- 0.5
			uy[[i]]$eaf.outcome <- 0.5		
		}
	}

	out$candidate_instruments <- subset(bind_rows(ux), select=c(
		SNP, beta.exposure, se.exposure, id.exposure, exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, pval.exposure
	))

	message("Removing known outliers from candidate instruments")

	out2 <- subset(out$search, sig)
	out$candidate_instruments <- group_by(out$candidate_instruments, id.exposure) %>%
		do({
			x <- .
			y <- subset(out2, id.outcome == x$id.exposure[1])
			x <- subset(x, !SNP %in% y$SNP)
			x
		})

	message("MR of candidates on outcome")

	out$candidate_outcome <- subset(bind_rows(uy), select=c(
		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
	))

	# out$candidate_outcome_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_outcome))
	out$candidate_outcome_dat <- harmonise_data(out$candidate_instruments, out$candidate_outcome)
	out$candidate_outcome_mr <- suppressMessages(mr(out$candidate_outcome_dat, method_list="mr_ivw"))


	message("MR of candidates on exposure")

	out$candidate_exposure <- subset(bind_rows(ux), select=c(
		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
	))

	# out$candidate_exposure_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_exposure))
	out$candidate_exposure_dat <- harmonise_data(out$candidate_instruments, out$candidate_exposure)
	out$candidate_exposure_mr <- suppressMessages(mr(out$candidate_exposure_dat, method_list="mr_ivw"))

	message("Finalising")
	out$G_annotation <- G_annotation
	out$simulation <- list()
	out$simulation$no_outlier_flag <- no_outlier_flag
	out$simulation$nid <- nid
	out$simulation$ngx <- ngx
	out$simulation$ngu1 <- ngu1
	out$simulation$ngu2 <- ngu2
	out$simulation$nu2 <- nu2
	out$simulation$ngu3 <- ngu3
	out$simulation$vgx <- vgx
	out$simulation$vgu1 <- vgu1
	out$simulation$vgu2 <- vgu2
	out$simulation$vgu3 <- vgu3
	out$simulation$bxy <- bxy
	out$simulation$bu1x <- bu1x
	out$simulation$bu1y <- bu1y
	out$simulation$bxu3 <- bxu3
	out$simulation$bu3y <- bu3y
	out$simulation$vgxu2 <- vgxu2
	out$simulation$vu2y <- vu2y
	out$simulation$mininum_instruments <- mininum_instruments
	out$simulation$instrument_threshold <- instrument_threshold
	out$simulation$outlier_threshold <- outlier_threshold
	out$simulation$outliers_known <- outliers_known
	out$simulation$directional_bias = directional_bias
	out$simulation$n_instruments <- sum(G_annotation$detected)
	out$simulation$n_valid_instruments <- sum(G_annotation$detected * G_annotation$valid)


return(out)
}









# #' Simulate data to test tryx
# #' 
# #' Choose a sample size, number of horizontal pleiotropy paths, and number of confounders. Create summary data for all necessary traits, then perform tryx scan
# #' 
# #' @param nid Sample size to simulate
# #' @param nu1 Number of traits which influence both X and Y
# #' @param nu2 Number of traits which influence Y
# #' @param bxu3 Effect of x on mediator
# #' @param bu3y Effect of mediator on y
# #' @param bxy Causal effect of X on Y
# #' @param outliers_known Assume that all outliers are detected (default = 'known'). But can also be 'detected', where we use heterogeneity based outlier detection, or 'all' where we try to adjust for every outlier regardless of heterogeneity
# #' @param directional_bias Is the outlier typically in one direction (i.e. unbalanced pleiotropy?)
# #' @param outlier_threshold Either 'nominal' or 'bonferroni'
# #' @param single_pleiotropy_trait = FALSE
# #' 
# #' @export
# #' @return Results from simulations and tryx scan
# tryx.simulate <- function(nid, nu1, nu2, bxu3=0, bu3y=0, bxy=3, outliers_known=c("known", "detected", "all")[1], directional_bias=FALSE, outlier_threshold = 'nominal', single_pleiotropy_trait = FALSE, debug=FALSE, ngu4=10, ngx=30, bu4x=0, bu4y=0)
# {
# 	# scenario 1 - confounder g -> u; u -> x; u -> y
# 	# scenario 2 - pleiotropy g -> u -> y
# 	# scenario 3 - mediator
# 	# scenario 4 - 

# 	out <- list()

# 	cpg <- require(simulateGP)
# 	if(!cpg)
# 	{
# 		stop("Please install the simulateGP package\ndevtools::install_github('explodecomputer/simulateGP')")
# 	}

# 	ngu1 <- 30
# 	ngu2 <- 30
# 	ngu3 <- 30
# 	nu3 <- 1
# 	nu4 <- 1
# 	nu <- nu1 + nu2 + nu3
# 	stopifnot(ngx >= nu)

# 	gx <- make_geno(nid, ngx, 0.5)
# 	out$true_outliers <- 1:nu

# 	# Effect sizes
# 	bgu1 <- lapply(1:nu1, function(x) runif(ngu1))
# 	bgu2 <- lapply(1:nu2, function(x) runif(ngu2))
# 	bgu3 <- runif(ngu3)
# 	bgx <- runif(ngx)

# 	bu1x <- out$bu1x <- rnorm(nu1)
# 	bu1y <- out$bu1y <- rnorm(nu1)
# 	bu2y <- out$bu2y <- rnorm(nu2)
# 	if(directional_bias)
# 	{
# 		bu1x <- out$bu1x <- abs(bu1x)
# 		bu1y <- out$bu1y <- abs(bu1y)
# 		bu2y <- out$bu2y <- abs(bu2y)
# 	}
# 	out$bxu3 <- bxu3
# 	bu3y <- out$bu3y <- bu3y
# 	bxy <- out$bxy <- bxy

# 	# create outlier traits
# 	gu1 <- list()
# 	u1 <- matrix(rnorm(nid * ngx), nid, ngx)
# 	if(nu1 > 0)
# 	{
# 		for(i in 1:nu1)
# 		{
# 			gu1[[i]] <- cbind(gx[,i], make_geno(nid, ngu1-1, 0.5))
# 			if(single_pleiotropy_trait)
# 			{
# 				u1[,1] <- u1[,1] + gu1[[i]] %*% bgu1[[i]]
# 			} else {
# 				u1[,i] <- gu1[[i]] %*% bgu1[[i]] + rnorm(nid)
# 			}
# 		}
# 	}

# 	gu2 <- list()
# 	u2 <- matrix(rnorm(nid * ngx), nid, ngx)
# 	if(nu2 > 0)
# 	{
# 		for(j in 1:nu2)
# 		{
# 			gu2[[j]] <- cbind(gx[,i+j], make_geno(nid, ngu1-1, 0.5))
# 			if(single_pleiotropy_trait)
# 			{
# 				u2[,1] <- u2[,1] + gu2[[j]] %*% bgu2[[j]]
# 			} else {
# 				u2[,j] <- gu2[[j]] %*% bgu2[[j]] + rnorm(nid)
# 			}
# 		}
# 	}

# 	u4 <- rnorm(nid)
# 	if(ngu4 > 0)
# 	{
# 		u4 <- u4 + gx[,1:ngu4] %*% runif(ngu4)
# 	}
# 	x <- drop(gx[,c((ngu4+1):ngx)] %*% bgx[c((ngu4+1):ngx)] + u4 %*% bu4x + u1 %*% bu1x + rnorm(nid))
	
# 	gu3 <- list(make_geno(nid, ngu3, 0.5))
# 	u3 <- matrix(gu3[[1]] %*% bgu3 + x * bxu3 + rnorm(nid), nid, 1)

# 	y <- drop(x * bxy + u1 %*% bu1y + u2 %*% bu2y + u3 * bu3y + u4 %*% bu4y + rnorm(nid))

# 	out$dat <- get_effs(x, y, gx, "X", "Y")
# 	out$dat$mr_keep <- TRUE
# 	no_outlier_flag <- FALSE

# 	if(outliers_known == "known")
# 	{
# 		message("Assuming knowledge of all outliers")
# 		outliers <- as.character(c(1:nu))
# 	} else if(outliers_known %in% c("detected", "all")) {
# 		message("Detecting outliers")
# 		radial <- RadialMR::ivw_radial(RadialMR::format_radial(out$dat$beta.exposure, out$dat$beta.outcome, out$dat$se.exposure, out$dat$se.outcome, out$dat$SNP), ifelse(outlier_threshold == "nominal", 0.05, 0.05/nrow(out$dat)))
# 		if(radial$outliers[1] == "No significant outliers")
# 		{
# 			message("No significant outliers detected")
# 			outliers <- as.character(c(1:nu))
# 			no_outlier_flag <- TRUE
# 		} else {
# 			outliers <- as.character(sort(radial$outliers$SNP))
# 			message("Number of outliers: ", length(outliers))
# 		}
# 	} else {
# 		stop("outliers_known must be known, detected or all")
# 	}
# 	# output$radialmr <- radial
# 	nout <- length(outliers)
# 	outliers <- subset(out$dat, SNP %in% outliers)$SNP
# 	out$outliers <- outliers
# 	if(outliers_known == "all")
# 	{
# 		message("Treating all variants as outliers")
# 		outliers <- 1:ngx
# 	}
# 	out$id_list <- c(paste0("U1_", 1:nu1), paste0("U2_", 1:nu2))

# 	# temp1 <- gwas(u1, gx[,1:2])

# 	temp <- list()
# 	for(i in 1:nu1)
# 	{
# 		temp[[i]] <- gwas(u1[,i], gx[,outliers, drop=FALSE])
# 		temp[[i]]$SNP <- outliers
# 		temp[[i]]$outcome <- paste0("U1_", i)
# 	}
# 	for(i in 1:nu2)
# 	{
# 		temp[[i + nu1]] <- gwas(u2[,i], gx[,outliers, drop=FALSE])
# 		temp[[i + nu1]]$SNP <- outliers
# 		temp[[i + nu1]]$outcome <- paste0("U2_", i)
# 	}
# 	for(i in 1:nu3)
# 	{	
# 		temp[[i + nu1 + nu2]] <- gwas(u3[,i], gx[,outliers, drop=FALSE])
# 		temp[[i + nu1 + nu2]]$SNP <- outliers
# 		temp[[i + nu1 + nu2]]$outcome <- paste0("U3_", i)
# 	}
# 	for(i in 1:nu4)
# 	{	
# 		temp[[i + nu1 + nu2 + nu3]] <- gwas(u4[,i], gx[,outliers, drop=FALSE])
# 		temp[[i + nu1 + nu2 + nu3]]$SNP <- outliers
# 		temp[[i + nu1 + nu2 + nu3]]$outcome <- paste0("U4_", i)
# 	}

# 	temp <- bind_rows(temp)

# 	out$search <- tibble::data_frame(SNP=temp$SNP, beta.outcome=temp$bhat, se.outcome=temp$se, pval.outcome=temp$pval, id.outcome=temp$outcome, outcome=temp$outcome, mr_keep=TRUE, effect_allele.outcome="A", other_allele.outcome="G", eaf.outcome=0.5)
# 	out$search$sig <- out$search$pval.outcome < 5e-8

# 	## Multivariable MR
# 	# Do mvmr with just the outlier detected traits

# 	traitlist <- as.character(unique(subset(out$search, sig)$id.outcome))
# 	G <- list()
# 	PHEN <- list(x)
# 	j <- 2
# 	if(any(grepl("U1", traitlist)))
# 	{
# 		u1traits <- gsub("U1_", "", traitlist[grepl("U1", traitlist)]) %>% as.numeric()
# 		message("U1 traits: ", paste(u1traits, collapse=","))
# 		G[[1]] <- lapply(u1traits, function(x) {
# 			return(gu1[[x]][,-1])
# 		}) %>% do.call(cbind, .)
# 		for(i in u1traits)
# 		{
# 			PHEN[[j]] <- u1[,i]
# 			j <- j + 1
# 		}
# 	}
# 	if(any(grepl("U2", traitlist)))
# 	{
# 		u2traits <- gsub("U2_", "", traitlist[grepl("U2", traitlist)]) %>% as.numeric()
# 		message("U2 traits: ", paste(u2traits, collapse=","))
# 		G[[2]] <- lapply(u2traits, function(x) {
# 			return(gu2[[x]][,-1])
# 		}) %>% do.call(cbind, .)
# 		for(i in u2traits)
# 		{
# 			PHEN[[j]] <- u2[,i]
# 			j <- j + 1
# 		}
# 	}
# 	if(any(grepl("U3", traitlist)))
# 	{
# 		u3traits <- gsub("U3_", "", traitlist[grepl("U3", traitlist)]) %>% as.numeric()
# 		message("U3 traits: ", paste(u3traits, collapse=","))
# 		G[[3]] <- lapply(u3traits, function(x) {
# 			return(gu3[[x]])
# 		}) %>% do.call(cbind, .)
# 		for(i in u3traits)
# 		{
# 			PHEN[[j]] <- u3[,i]
# 			j <- j + 1
# 		}
# 	}
# 	if(any(grepl("U4", traitlist)))
# 	{
# 		u4traits <- gsub("U4_", "", traitlist[grepl("U4", traitlist)]) %>% as.numeric()
# 		message("U4 traits: ", paste(u4traits, collapse=","))
# 		G[[4]] <- lapply(u4traits, function(x) {
# 			return(gx)
# 		}) %>% do.call(cbind, .)
# 		for(i in u4traits)
# 		{
# 			PHEN[[j]] <- u4[,i]
# 			j <- j + 1
# 		}
# 	}
# 	G[[5]] <- gx
# 	G <- do.call(cbind, G)
# 	traitlist <- c("X", traitlist)
# 	names(PHEN) <- traitlist

# 	out$mvres <- try({
# 		mvdat <- make_mvdat(PHEN, y, G)
# 		mvres <- mv_multiple(mvdat)
# 		mvres$result$exposure <- traitlist
# 		mvres
# 	})
# 	## get effects
# 	u1x <- list()
# 	u1y <- list()
# 	for(i in 1:nu1)
# 	{
# 		u1x[[i]] <- get_effs(u1[,i], x, gu1[[i]], paste0("U1_", i), "X")
# 		u1y[[i]] <- get_effs(u1[,i], y, gu1[[i]], paste0("U1_", i), "Y")
# 		u1x[[i]]$SNP <- u1y[[i]]$SNP <- c(i, paste0("u1_",1:(ngu1-1)))
# 		u1x[[i]]$effect_allele.exposure <- u1y[[i]]$effect_allele.exposure <- "A"
# 		u1x[[i]]$other_allele.exposure <- u1y[[i]]$other_allele.exposure <- "G"
# 		u1x[[i]]$effect_allele.outcome <- u1y[[i]]$effect_allele.outcome <- "A"
# 		u1x[[i]]$other_allele.outcome <- u1y[[i]]$other_allele.outcome <- "G"
# 		u1x[[i]]$eaf.exposure <- u1y[[i]]$eaf.exposure <- 0.5
# 		u1x[[i]]$eaf.outcome <- u1y[[i]]$eaf.outcome <- 0.5
# 	}

# 	u2x <- list()
# 	u2y <- list()
# 	for(i in 1:nu2)
# 	{
# 		u2x[[i]] <- get_effs(u2[,i], x, gu2[[i]], paste0("U2_", i), "X")
# 		u2y[[i]] <- get_effs(u2[,i], y, gu2[[i]], paste0("U2_", i), "Y")
# 		u2x[[i]]$SNP <- u2y[[i]]$SNP <- c(i+nu1, paste0("u2_",1:(ngu2-1)))
# 		u2x[[i]]$effect_allele.exposure <- u2y[[i]]$effect_allele.exposure <- "A"
# 		u2x[[i]]$other_allele.exposure <- u2y[[i]]$other_allele.exposure <- "G"
# 		u2x[[i]]$effect_allele.outcome <- u2y[[i]]$effect_allele.outcome <- "A"
# 		u2x[[i]]$other_allele.outcome <- u2y[[i]]$other_allele.outcome <- "G"
# 		u2x[[i]]$eaf.exposure <- u2y[[i]]$eaf.exposure <- 0.5
# 		u2x[[i]]$eaf.outcome <- u2y[[i]]$eaf.outcome <- 0.5
# 	}

# 	u3x <- list()
# 	u3y <- list()
# 	for(i in 1:nu3)
# 	{
# 		u3x[[i]] <- get_effs(u3[,i], x, gu3[[i]], paste0("U3_", i), "X")
# 		u3y[[i]] <- get_effs(u3[,i], y, gu3[[i]], paste0("U3_", i), "Y")
# 		u3x[[i]]$SNP <- u3y[[i]]$SNP <- c(i+nu1, paste0("u3_",1:(ngu3-1)))
# 		u3x[[i]]$effect_allele.exposure <- u3y[[i]]$effect_allele.exposure <- "A"
# 		u3x[[i]]$other_allele.exposure <- u3y[[i]]$other_allele.exposure <- "G"
# 		u3x[[i]]$effect_allele.outcome <- u3y[[i]]$effect_allele.outcome <- "A"
# 		u3x[[i]]$other_allele.outcome <- u3y[[i]]$other_allele.outcome <- "G"
# 		u3x[[i]]$eaf.exposure <- u3y[[i]]$eaf.exposure <- 0.5
# 		u3x[[i]]$eaf.outcome <- u3y[[i]]$eaf.outcome <- 0.5
# 	}

# 	u4x <- list()
# 	u4y <- list()
# 	for(i in 1:nu4)
# 	{
# 		u4x[[i]] <- get_effs(u4[,i], x, gx[,1:ngu4], paste0("U4_", i), "X")
# 		u4y[[i]] <- get_effs(u4[,i], y, gx[,1:ngu4], paste0("U4_", i), "Y")
# 		u4x[[i]]$SNP <- u4y[[i]]$SNP <- c(i+nu1, paste0("u4_",1:(ngu4-1)))
# 		u4x[[i]]$effect_allele.exposure <- u4y[[i]]$effect_allele.exposure <- "A"
# 		u4x[[i]]$other_allele.exposure <- u4y[[i]]$other_allele.exposure <- "G"
# 		u4x[[i]]$effect_allele.outcome <- u4y[[i]]$effect_allele.outcome <- "A"
# 		u4x[[i]]$other_allele.outcome <- u4y[[i]]$other_allele.outcome <- "G"
# 		u4x[[i]]$eaf.exposure <- u4y[[i]]$eaf.exposure <- 0.5
# 		u4x[[i]]$eaf.outcome <- u4y[[i]]$eaf.outcome <- 0.5
# 	}


# 	out$candidate_instruments <- subset(rbind(bind_rows(u1x), bind_rows(u2x), bind_rows(u3x)), select=c(
# 		SNP, beta.exposure, se.exposure, id.exposure, exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, pval.exposure
# 	))
# 	out2 <- subset(out$search, sig)
# 	out$candidate_instruments <- group_by(out$candidate_instruments, id.exposure) %>%
# 		do({
# 			x <- .
# 			y <- subset(out2, id.outcome == x$id.exposure[1])
# 			x <- subset(x, !SNP %in% y$SNP)
# 			x
# 		})

# 	out$candidate_outcome <- subset(rbind(bind_rows(u1y), bind_rows(u2y), bind_rows(u3y)), select=c(
# 		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
# 	))

# 	# out$candidate_outcome_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_outcome))
# 	out$candidate_outcome_dat <- bind_rows(bind_rows(u1y), bind_rows(u2y), bind_rows(u3y))
# 	out$candidate_outcome_mr <- suppressMessages(mr(out$candidate_outcome_dat, method_list="mr_ivw"))


# 	out$candidate_exposure <- subset(rbind(bind_rows(u1x), bind_rows(u2x), bind_rows(u3x)), select=c(
# 		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
# 	))

# 	# out$candidate_exposure_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_exposure))
# 	out$candidate_exposure_dat <- bind_rows(bind_rows(u1x), bind_rows(u2x), bind_rows(u3x))
# 	out$candidate_exposure_mr <- suppressMessages(mr(out$candidate_exposure_dat, method_list="mr_ivw"))


# 	out$simulation <- list()
# 	out$simulation$no_outlier_flag <- no_outlier_flag
# 	out$simulation$nu1 <- nu1
# 	out$simulation$nu2 <- nu2
# 	out$simulation$nu3 <- nu3
# 	out$simulation$bxu3 <- bxu3
# 	out$simulation$bu3y <- bu3y
# 	out$simulation$outliers_known <- outliers_known
# 	out$simulation$bxy <- bxy

# 	if(debug)
# 	{
# 		out$phen <- list(
# 			y=y, x=x, u1=u1, u2=u2, u3=u3, u4=u4
# 		)
# 	}
# 	return(out)
# }

