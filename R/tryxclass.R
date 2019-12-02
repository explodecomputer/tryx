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
	},
  # Find associations with outliers
  set_candidate_traits = function(id_list=NULL) {
  	if(is.null(id_list))
  	{
        ao <- suppressMessages(TwoSampleMR::available_outcomes())
        ids <- subset(ao, priority == 1 & nsnp > 500000 & sample_size > 5000) %>%
               arrange(desc(sample_size)) %>%
               filter(!duplicated(trait), mr == 1) %>%
               filter(!author %in% c("Shin", "Kettunen", "Roederer")) %>%
               filter(!grepl("ukb-b", id)) %>%
               filter(! id %in% c(dat$id.exposure[1], dat$id.outcome[1]))
        id_list <- ids$id
        message("Using default list of ", nrow(ids), " traits")
   	 }
  	self$output$idlist <- id_list
    	invisible(self$output$idlist)
 	},
  #Search for candidate traits associated with outliers
  scan_candidate_traits = function(dat=self$output$dat, search_correction="none", search_threshold=ifelse(search_correction=="none", 5e-8, 0.05), use_proxies=FALSE) {
	    stopifnot(search_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
 	    stopifnot(search_threshold > 0 & search_threshold < 1)
 	    self$output$search <- extract_outcome_data(self$output$outliers, self$output$id_list, proxies=use_proxies)
 	    self$output$search$padj <- p.adjust(self$output$search$pval.outcome, search_correction)
    	    out2 <- subset(self$output$search, padj < search_threshold)
   	 if(nrow(out2) == 0)
         {
     	  message("Outliers do not associate with any other traits. Try relaxing the search_threshold")
   	  return(self$output)
     	  }
    	    message("Found ", length(unique(out2$id.outcome)), " candidate traits associated with outliers at p < ", search_threshold)
            invisible(self$output$search)
    	    self$output$sig_search <- out2
    	    invisible(self$output$sig_search)
  	},
  #Obtain instruments for the candidate traits
  candidate_instruments = function(candidate_instruments = NULL, include_outliers = FALSE) {
  	    message("Finding instruments for candidate traits")
  	    if(is.null(candidate_instruments))
  	    {
         	candidate_instruments <- extract_instruments(unique(self$output$sig_search$id.outcome))
         	self$output$candidate_instruments <- candidate_instruments
         	if(nrow(candidate_instruments) == 0)
         	{
            	message("No instruments available for the candidate traits")
            	return(self$output)
         	}
                if(!include_outliers)
       	        {
           	message("Removing outlier SNPs from candidate trait instrument lists")
            	candidate_instruments <- group_by(candidate_instruments, id.exposure) %>%
                                      do({
                                          z <- .
                                          y <- subset(self$output$sig_search, id.outcome == z$id.exposure[1])
                                          z <- subset(z, !SNP %in% y$SNP)
                                          z
                                        })
                 }
  	       }
  	    invisible(self$output$candidate_instruments)
  	    message(length(unique(self$output$candidate_instruments$id.exposure)), " traits with at least one instrument")
  	 }
))
  
