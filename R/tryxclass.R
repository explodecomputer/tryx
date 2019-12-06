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
    self$output$id_list <- id_list
    invisible(self$output$id_list)
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
  message(length(unique(self$output$candidate_instruments$id.exposure)), " traits with at least one instrument")
  invisible(self$output$candidate_instruments)
  },
  #extract instrument for candidate trait instruments for the original outcome
  outcome_instruments = function(candidate_outcome = NULL, dat = self$output$dat, use_proxies=FALSE) {
    message("Looking up candidate trait instruments for ", dat$outcome[1])
    candidate_outcome <- extract_outcome_data(unique(self$output$candidate_instruments$SNP), dat$id.outcome[1], proxies=use_proxies)
    self$output$candidate_outcome <- candidate_outcome
    if(is.null(candidate_outcome))
    {
        message("None of the candidate trait instruments available for ", dat$outcome[1])
        return(self$output)
    }
    message(nrow(candidate_outcome), " instruments extracted for ", dat$outcome[1])
    invisible(self$output$candidate_outcome)
  },
  #make a dataset for the candidate traits and the original outcome
  candidate_outcome_dat = function(dat = self$output$dat){
    candidate_outcome_dat <- suppressMessages(harmonise_data(self$output$candidate_instruments, self$output$candidate_outcome))
    candidate_outcome_dat <- subset(candidate_outcome_dat, mr_keep)
    if(nrow(candidate_outcome_dat) == 0)
    {
      message("None of the candidate trait instruments available for ", dat$outcome[1], " after harmonising")
      return(self$output)
    }
    self$output$candidate_outcome_dat <- candidate_outcome_dat
    invisible(self$output$candidate_outcome_dat)
  },
  #extract instrument for candidate trait instruments for the original exposure
  exposure_instruments = function(candidate_exposure = NULL, dat = self$output$dat, use_proxies=FALSE) {
    message("Looking up candidate trait instruments for ", dat$exposure[1])
    candidate_exposure <- extract_outcome_data(unique(self$output$candidate_instruments$SNP), dat$id.exposure[1], proxies=use_proxies)
    self$output$candidate_exposure <- candidate_exposure
    if(is.null(candidate_exposure))
    {
      message("None of the candidate trait instruments available for ", dat$exposure[1])
      return(self$output)
    }
    message(nrow(self$output$candidate_exposure), " instruments extracted for ", dat$exposure[1])
    invisible(self$output$candidate_exposure)
  },
  #make a dataset for the candidate traits and the original exposure
  candidate_exposure_dat = function(dat = self$output$dat){
    candidate_exposure_dat <- suppressMessages(harmonise_data(self$output$candidate_instruments, self$output$candidate_exposure))
    candidate_exposure_dat <- subset(candidate_exposure_dat, mr_keep)
    self$output$candidate_exposure_dat <- candidate_exposure_dat
    if(nrow(candidate_exposure_dat) == 0)
    {
      message("None of the candidate trait instruments available for ", dat$exposure[1], " after harmonising")
      return(self$output)
    }
  },
  #extract instrument for the original exposure for candidate traits
  exposure_candidate_instruments = function(exposure_candidate = NULL, dat = self$output$dat, use_proxies = FALSE, include_outliers = FALSE) {
    message("Looking up exposure instruments for ", length(unique(self$output$sig_search$id.outcome)), " candidate traits")
    exposure_candidate <- extract_outcome_data(unique(dat$SNP), unique(self$output$sig_search$id.outcome), proxies=use_proxies)
    self$output$exposure_candidate <- exposure_candidate
    if(is.null(exposure_candidate))
    {
      message("None of the candidate trait instruments available for ", dat$exposure[1])
      return(self$output)
    }
    message(nrow(self$output$candidate_exposure), " instruments extracted")
    temp <- subset(dat, select=c(SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, id.exposure, exposure))
    if(!include_outliers)
    {
      message("Removing outlier SNPs from candidate trait outcome lists")
      n1 <- nrow(self$output$exposure_candidate)
      self$output$exposure_candidate <- group_by(self$output$exposure_candidate, id.outcome) %>%
      do({
        z <- .
        y <- subset(self$output$sig_search, id.outcome == z$id.outcome[1])
        z <- subset(z, !SNP %in% y$SNP)
        z
      })
      message("Removed ", n1 - nrow(self$output$exposure_candidate), " outlier SNPs")
    }
    invisible(self$output$exposure_candidate)
    self$output$temp <- temp
    invisible(self$output$temp)
  },
  #make a dataset for the original exposure and the candidate traits
  exposure_candidate_dat = function(dat = self$output$dat){
     exposure_candidate_dat <- suppressMessages(harmonise_data(self$output$temp, self$output$exposure_candidate))
     exposure_candidate_dat <- subset(exposure_candidate_dat, mr_keep)
     self$output$exposure_candidate_dat <- exposure_candidate_dat
     if(nrow(exposure_candidate_dat) == 0)
     {
        message("None of the candidate traits have the ", dat$exposure[1], " instruments after harmonising")
        return(self$output)
     }
     invisible(self$output$exposure_candidate_dat)
  },
  #Performing MR of 1) candidate traits-outcome 2) candidate traits-exposure 3) exposure-candidate traits
  perform_mr = function(dat = self$output$dat, mr_method="mr_ivw") {
    message("Performing MR of ", length(unique(self$output$candidate_outcome_dat$id.exposure)), " candidate traits against ", dat$outcome[1])
    if(mr_method == "strategy1")
    {
      temp <- tryx::strategy1(self$output$candidate_outcome_dat)
      self$output$candidate_outcome_mr <- temp$res
      self$output$candidate_outcome_mr_full <- temp
    } else {
      self$output$candidate_outcome_mr <- suppressMessages(mr(self$output$candidate_outcome_dat, method_list=c("mr_wald_ratio", mr_method)))
    }
    invisible(self$output$candidate_outcome_mr)
    invisible(self$output$candidate_outcome_mr_full)

    message("Performing MR of ", length(unique(self$output$candidate_exposure_dat$id.exposure)), " candidate traits against ", dat$exposure[1])
    if(mr_method == "strategy1")
    {
      temp <- tryx::strategy1(self$output$candidate_exposure_dat)
      self$output$candidate_exposure_mr <- temp$res
      self$output$candidate_exposure_mr_full <- temp
    } else {
      self$output$candidate_exposure_mr <- suppressMessages(mr(self$output$candidate_exposure_dat, method_list=c("mr_wald_ratio", mr_method)))
    }
    invisible(self$output$candidate_exposure_mr)
    invisible(self$output$candidate_exposure_mr_full)
    
    message("Performing MR of ", dat$exposure[1], " against ", length(unique(self$output$exposure_candidate_dat$id.outcome)), " candidate traits")
    if(mr_method == "strategy1")
    {
      temp <- tryx::strategy1(self$output$exposure_candidate_dat)
      self$output$exposure_candidate_mr <- temp$res
      self$output$exposure_candidate_mr_full <- temp
    } else {
      self$output$exposure_candidate_mr <- suppressMessages(mr(self$output$exposure_candidate_dat, method_list=c("mr_wald_ratio", mr_method)))
    }
    invisible(self$output$exposure_candidate_mr)
    invisible(self$output$exposure_candidate_mr_full)
  },

  # GIB BROKE THIS

  scan = function(mr_method="mr_ivw") {
    self$
    self$
    self$perform_mr(mr_method=mr_method)
  }


))
 
