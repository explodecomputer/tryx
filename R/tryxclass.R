Tryx <- R6::R6Class("Tryx", list(
	output = list(),
	input = function(dat) {
		if(length(unique(dat$id.exposure)) > 1 | length(unique(dat$id.outcome)) > 1)
		{
			message("Warning! Multiple exposure/outcome combinations found")
			message("Only using first exposure / outcome combination")
		}
		dat <- subset(dat, id.exposure == id.exposure[1] & id.outcome == id.outcome[1] & mr_keep)
		self$output$dat <- dat
	    invisible(self)
	},
	test = function(dat = self$output$dat) {
		print(dat)
		invisible(NULL)
	},
get_outliers = function(dat=self$output$dat, outliers="RadialMR", outlier_correction="none", outlier_threshold=ifelse(outlier_correction=="none", 0.05/nrow(dat), 0.05)) {
	stopifnot(outlier_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
	stopifnot(outlier_threshold > 0 & outlier_threshold < 1)
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
		self$output$radialmr <- radial
		if(length(outliers) == 0)
		{
			message("No outliers found. Exiting")
			return(self$output)
		}
	} else {
		nout <- length(outliers)
		outliers <- subset(dat, SNP %in% outliers)$SNP
		message(length(outliers), " of ", nout, " of specified outliers present in dat")
		if(length(outliers) == 0)
		{
			message("No outliers found. Exiting")
			return(self$output)
		}
	}
	self$output$outliers <- outliers
	invisible(self$output$outliers)
}

))
