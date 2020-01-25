#' Class for MR-TRYX analysis
#'
#' @description
#' A simple wrapper function.
#' Using a summary set, find outliers in the MR analysis between the pair of trais.
#' Find other 'candidate traits' associated with those outliers.
#' Perform MR of each of those candidate traits with the original exposure and outcome.
#' @export
Tryx <- R6::R6Class("Tryx", list(
  output = list(),
  
  ##########################################################################################################################################################################

#' @description
#' Create a new dataset and initialise an R interface
#' @param dat Dataset from TwoSampleMR::harmonise_data
#' 

  initialize = function(dat) {
    if(length(unique(dat$id.exposure)) > 1 | length(unique(dat$id.outcome)) > 1)
    {
      message("Warning! Multiple exposure/outcome combinations found")
      message("Only using first exposure / outcome combination")
    }
    dat <- subset(dat, id.exposure == id.exposure[1] & id.outcome == id.outcome[1] & mr_keep)
    self$output$dat <- dat
    invisible(self)
  },

  print = function(...) {
    cat("Tryx analysis of ", self$output$dat$exposure[1], " against ", self$output$dat$outcome[1], "\n")
    cat("Status:\n")
    null <- sapply(names(self$output), function(x) cat("  ", x, "\n"))
    invisible(self)
  },
  
  ##########################################################################################################################################################################
# Get a list of outliers used 
#' @description 
#' Detect outliers in exposure-outcome dataset.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.
#' 
#' @param outliers Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers.
#' 
#' @param outlier_correction Defualt = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#'  
#' @param outlier_threshold If outlier_correction = "none" then the p-value threshold for detecting outliers is by default 0.05.
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
  
  ##########################################################################################################################################################################
  # Find associations with outliers
#' @description 
#' Set a list of GWAS IDs used.
#' 
#' @param id_list The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes().
#' 
  set_candidate_traits = function(id_list=NULL) {
    if(is.null(id_list))
    {
      ao <- suppressMessages(TwoSampleMR::available_outcomes())
      ids <- subset(ao) %>% 
             arrange(desc(sample_size)) %>%
             filter(!duplicated(trait), mr == 1) %>%
             filter(!grepl("ukb-a", id)) %>%
             filter(! id %in% c(dat$id.exposure[1], dat$id.outcome[1]))
      id_list <- ids$id
      message("Using default list of ", nrow(ids), " traits")
    }
    self$output$id_list <- id_list
    invisible(self$output$id_list)
  },
  
  ##########################################################################################################################################################################
#' @description 
#' Search for candidate traits associated with outliers.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. 
#' 
#' @param search_correction Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none.
#' 
#' @param search_threshold If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.
#' 
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
  scan = function(dat=self$output$dat, search_correction="none", search_threshold=ifelse(search_correction=="none", 5e-8, 0.05), use_proxies=FALSE) {
    stopifnot(search_correction %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    stopifnot(search_threshold > 0 & search_threshold < 1)
    search <- extract_outcome_data(self$output$outliers, self$output$id_list, proxies=use_proxies)
    padj <- p.adjust(search$pval.outcome, search_correction)
    search$sig <- padj < search_threshold
    out2 <- subset(search, sig)
    if(nrow(out2) == 0)
    {
      message("Outliers do not associate with any other traits. Try relaxing the search_threshold")
      return(self$output)
    }
    message("Found ", length(unique(out2$id.outcome)), " candidate traits associated with outliers at p < ", search_threshold)
    self$output$search <- search
    invisible(self$output$search)
    self$output$sig_search <- out2
    invisible(self$output$sig_search)
  },
  
##########################################################################################################################################################################
#Instruments for candidate traits 
#' @description
#' Obtain instruments for the candidate traits.
#' 
#' @param candidate_instruments Instruments for candidate traits.
#' 
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE. 
  candidate_instruments = function(candidate_instruments = NULL, include_outliers = FALSE) {
  message("Finding instruments for candidate traits")
  if(is.null(candidate_instruments))
  {
         candidate_instruments <- extract_instruments(unique(self$output$sig_search$id.outcome))
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
  message(length(unique(candidate_instruments$id.exposure)), " traits with at least one instrument")
  self$output$candidate_instruments <- candidate_instruments
  invisible(self$output$candidate_instruments)
  },
  
##########################################################################################################################################################################
#Instruments for the outcome
#' @description 
#' Extract instrument for candidate trait instruments for the original outcome.
#' 
#' @param candidate_outcome Extracted instrument SNPs from outcome.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data.
#' 
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
  outcome_instruments = function(candidate_outcome = NULL, dat = self$output$dat, use_proxies=FALSE) {
    message("Looking up candidate trait instruments for ", dat$outcome[1])
    candidate_outcome <- extract_outcome_data(unique(self$output$candidate_instruments$SNP), dat$id.outcome[1], proxies=use_proxies)
    if(is.null(candidate_outcome))
    {
        message("None of the candidate trait instruments available for ", dat$outcome[1])
        return(self$output)
    }
    message(nrow(candidate_outcome), " instruments extracted for ", dat$outcome[1])
    self$output$candidate_outcome <- candidate_outcome
    invisible(self$output$candidate_outcome)
  },

########################################################################################################################################################################## 
#Instruments for the exposure 
#' @description
#' Extract instrument for candidate trait instruments for the original exposure.
#' 
#' @param candidate_exposure Extracted instrument SNPs from exposure.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data.
#' 
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
  exposure_instruments = function(candidate_exposure = NULL, dat = self$output$dat, use_proxies=FALSE) {
    message("Looking up candidate trait instruments for ", dat$exposure[1])
    candidate_exposure <- extract_outcome_data(unique(self$output$candidate_instruments$SNP), dat$id.exposure[1], proxies=use_proxies)
    if(is.null(candidate_exposure))
    {
      message("None of the candidate trait instruments available for ", dat$exposure[1])
      return(self$output)
    }
    message(nrow(self$output$candidate_exposure), " instruments extracted for ", dat$exposure[1])
    self$output$candidate_exposure <- candidate_exposure
    invisible(self$output$candidate_exposure)
  },

########################################################################################################################################################################## 
#Extracted instrument SNPs from exposure 
#' @description 
#' Extract instrument for the original exposure for the candidate traits.
#' 
#' @param exposure_candidate Extracted instrument SNPs from exposure.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data.
#' 
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
#' 
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
  exposure_candidate_instruments = function(exposure_candidate = NULL, dat = self$output$dat, use_proxies = FALSE, include_outliers = FALSE) {
    message("Looking up exposure instruments for ", length(unique(self$output$sig_search$id.outcome)), " candidate traits")
    exposure_candidate <- extract_outcome_data(unique(dat$SNP), unique(self$output$sig_search$id.outcome), proxies=use_proxies)
    if(is.null(exposure_candidate))
    {
      message("None of the candidate trait instruments available for ", dat$exposure[1])
      return(self$output)
    }
    self$output$exposure_candidate <- exposure_candidate
    invisible(self$output$exposure_candidate)
    message(nrow(self$output$exposure_candidate), " instruments extracted")
    
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
    self$output$temp <- temp
    invisible(self$output$temp)
  },
 
##########################################################################################################################################################################  
  # All in one - Extraction
#' @description 
#' Extract instruments for MR analyses.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data.
#'
#' @param candidate_instruments Instruments for candidate traits.
#'
#' @param candidate_outcome Extracted instrument SNPs from outcome.
#'
#' @param candidate_exposure Extracted instrument SNPs from exposure.
#'
#' @param exposure_candidate Extracted instrument SNPs from exposure.
#'
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE. 
#' 
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
  extractions = function(dat = self$output$dat, candidate_instruments = NULL, candidate_outcome = NULL, candidate_exposure = NULL, exposure_candidate = NULL, include_outliers = FALSE, use_proxies=FALSE){
    x$candidate_instruments(candidate_instruments = NULL, include_outliers = FALSE)
    x$outcome_instruments(candidate_outcome = NULL, dat = self$output$dat, use_proxies=FALSE)
    x$exposure_instruments(candidate_exposure = NULL, dat = self$output$dat, use_proxies=FALSE)
    x$exposure_candidate_instruments(exposure_candidate = NULL, dat = self$output$dat, use_proxies = FALSE, include_outliers = FALSE)
    invisible(self)
  },
  
  ##########################################################################################################################################################################
#Harmonised candidate traits - outcome dataset 
#' @description 
#' Make a dataset for the candidate traits and the original outcome.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. 
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
  
##########################################################################################################################################################################
#Harmonised candidate traits - exposure dataset 
#' @description 
#' Make a dataset for the candidate traits and the original exposure.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. 
  candidate_exposure_dat = function(dat = self$output$dat){
    candidate_exposure_dat <- suppressMessages(harmonise_data(self$output$candidate_instruments, self$output$candidate_exposure))
    candidate_exposure_dat <- subset(candidate_exposure_dat, mr_keep)
    if(nrow(candidate_exposure_dat) == 0)
    {
      message("None of the candidate trait instruments available for ", dat$exposure[1], " after harmonising")
      return(self$output)
    }
    self$output$candidate_exposure_dat <- candidate_exposure_dat
    invisible(self$output$candidate_exposure_dat)
  },

##########################################################################################################################################################################
#Harmonised exposure - candidate traits dataset 
#' @description 
#' Make a dataset for the original exposure and the candidate traits.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. 
  exposure_candidate_dat = function(dat = self$output$dat){
     exposure_candidate_dat <- suppressMessages(harmonise_data(self$output$temp, self$output$exposure_candidate))
     exposure_candidate_dat <- subset(exposure_candidate_dat, mr_keep)
     if(nrow(exposure_candidate_dat) == 0)
     {
        message("None of the candidate traits have the ", dat$exposure[1], " instruments after harmonising")
        return(self$output)
     }
     self$output$exposure_candidate_dat <- exposure_candidate_dat
     invisible(self$output$exposure_candidate_dat)
  },

##########################################################################################################################################################################
# All in one - data harmonisation
#' @description 
#' Harmonised exposure - outcome dataset. 
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. 
  harmonise = function(dat = self$output$dat) {
    x$candidate_outcome_dat(dat = self$output$dat)
    x$candidate_exposure_dat(dat = self$output$dat)
    x$exposure_candidate_dat(dat = self$output$dat)
    invisible(self)
  },
  
  ##########################################################################################################################################################################
  #MR analysis of candidates against exposure 
#' @description
#' Perform MR anlayses of 1) candidate traits-outcome 2) candidate traits-exposure 3) exposure-candidate traits.
#' 
#' @param dat Output from TwoSampleMR::harmonise_data. 
#' 
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is mr_ivw. Can also provide basic MR methods e.g. mr_weighted_mode, mr_weighted_median etc. Also possible to use "strategy1" which performs IVW in the first instance, but then weighted mode for associations with high heterogeneity.
  mr = function(dat = self$output$dat, mr_method="mr_ivw") {
    message("Performing MR of ", length(unique(self$output$candidate_outcome_dat$id.exposure)), " candidate traits against ", dat$outcome[1])
    if(mr_method == "strategy1")
    {
      temp <- x$strategy1(self$output$candidate_outcome_dat)
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
      temp <- x$strategy1(self$output$candidate_exposure_dat)
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
  
  ##########################################################################################################################################################################
  #All-in-one: Do everything at once
#' @description 
#' All-in-one: 1) Detect outlier 2) Search candidate traits 3) Perform MR of candidate traits and the outcome / exposure. 
#' 
#' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.
#' 
#' @param outliers Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers.
#' 
#' @param outlier_correction Defualt = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#' 
#' @param outlier_threshold If outlier_correction = "none" then the p-value threshold for detecting outliers is by default 0.05.
#' 
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed.
#' 
#' @param search_correction Default = "none", but can select from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#' 
#' @param search_threshold If search_correction = "none" then the p-value threshold for detecting an association between an outlier and a candidate trait is by default 5e-8. Otherwise it is 0.05.
#' 
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#' 
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is mr_ivw. Can also provide basic MR methods e.g. mr_weighted_mode, mr_weighted_median etc. Also possible to use "strategy1" which performs IVW in the first instance, but then weighted mode for associations with high heterogeneity.
  mrtryx = function(dat=self$output$dat, outliers="RadialMR", outlier_correction="none", outlier_threshold=ifelse(outlier_correction=="none", 0.05/nrow(dat), 0.05), use_proxies=FALSE, search_correction="none", search_threshold=ifelse(search_correction=="none", 5e-8, 0.05), include_outliers=FALSE, mr_method="mr_ivw") {
    x$get_outliers(dat=self$output$dat, outliers="RadialMR", outlier_correction="none", outlier_threshold=ifelse(outlier_correction=="none", 0.05/nrow(dat), 0.05))
    x$set_candidate_traits(id_list=NULL)
    x$scan(dat=self$output$dat, search_correction="none", search_threshold=ifelse(search_correction=="none", 5e-8, 0.05), use_proxies=FALSE)
    x$candidate_instruments(candidate_instruments = NULL, include_outliers = FALSE)
    x$outcome_instruments(candidate_outcome = NULL, dat = self$output$dat, use_proxies=FALSE)
    x$exposure_instruments(candidate_exposure = NULL, dat = self$output$dat, use_proxies=FALSE)
    x$exposure_candidate_instruments(exposure_candidate = NULL, dat = self$output$dat, use_proxies = FALSE, include_outliers = FALSE)
    x$candidate_outcome_dat(dat = self$output$dat)
    x$candidate_exposure_dat(dat = self$output$dat)
    x$exposure_candidate_dat(dat = self$output$dat)
    x$mr(dat = self$output$dat, mr_method="mr_ivw")
    invisible(self)
  },
  
  ##########################################################################################################################################################################
  #tryx-sig
#' @description
#' Identify putatively significant associations in the outlier scan.
#'  
#' @param mr_threshold_method This is the argument to be passed to \code{p.adjust}. Default is "fdr". If no p-value adjustment is to be applied then specify "unadjusted".
#' 
#' @param mr_threshold Threshold to declare significance  
  tryx.sig = function(mr_threshold_method = "fdr", mr_threshold = 0.05){

    stopifnot("candidate_outcome_mr" %in% names(self$output))
    stopifnot("candidate_exposure_mr" %in% names(self$output))
    
    # Use threshold to retain causal relationships
    stopifnot(length(mr_threshold_method) == 1)
    if(mr_threshold_method != "unadjusted")
    {
      message("Adjusting p-value")
      self$output$candidate_outcome_mr$pval_adj <- p.adjust(self$output$candidate_outcome_mr$pval, mr_threshold_method)
      self$output$candidate_outcome_mr$sig <- self$output$candidate_outcome_mr$pval_adj < mr_threshold
      self$output$candidate_exposure_mr$pval_adj <- p.adjust(self$output$candidate_exposure_mr$pval, mr_threshold_method)
      self$output$candidate_exposure_mr$sig <- self$output$candidate_exposure_mr$pval_adj < mr_threshold
      self$output$exposure_candidate_mr$pval_adj <- p.adjust(self$output$exposure_candidate_mr$pval, mr_threshold_method)
      self$output$exposure_candidate_mr$sig <- self$output$exposure_candidate_mr$pval_adj < mr_threshold
    } else {
      self$output$candidate_outcome_mr$sig <- self$output$candidate_outcome_mr$pval < mr_threshold
      self$output$candidate_exposure_mr$sig <- self$output$candidate_exposure_mr$pval < mr_threshold
      self$output$exposure_candidate_mr$sig <- self$output$exposure_candidate_mr$pval < mr_threshold
    }

    message("* * * *")
    message("Number of candidate - outcome associations: ", sum(self$output$candidate_outcome_mr$sig))
    message("* * * *")
    message(paste(subset(self$output$candidate_outcome_mr, sig)$exposure, collapse="\n"))
    message("")
    message("* * * *")
    message("Number of candidate - exposure associations: ", sum(self$output$candidate_exposure_mr$sig))
    message("* * * *")
    message(paste(subset(self$output$candidate_exposure_mr, sig)$exposure, collapse="\n"))
    message("")
    message("* * * *")
    message("Number of exposure - candidate associations: ", sum(self$output$exposure_candidate_mr$sig))
    message("* * * *")
    message(paste(subset(self$output$candidate_exposure_mr, sig)$exposure, collapse="\n"))
    #return(self$output)
    invisible(self)
  },


############################################################################################################################
# Functions for the adjustment method
############################################################################################################################
  #' @description
  #' Outlier adjustment estimation - How much of the heterogeneity due to the outlier can be explained by alternative pathways?
  #' 
  #' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.
  #' 
  #' @param tryxscan Output from \code{x$mrtryx()}
  #' 
  #' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
   adjustment = function(tryxscan=self$output, id_remove=NULL) {
    if(!any(tryxscan$search$sig))
    {
      return(NULL)
    }
    
    l <- list()
    sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
    sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)
    sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)
    
    dat <- tryxscan$dat
    dat$qi <- private$cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
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
          b <- private$bootstrap_path(a$beta.exposure, a$se.exposure, sig$beta.outcome[i], sig$se.outcome[i], sige$b[sige$id.exposure == sig$id.outcome[i]], sige$se[sige$id.exposure == sig$id.outcome[i]])
          a$adj.beta.exposure <- b[1]
          a$adj.se.exposure <- b[2]
          b <- private$bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
          a$adj.beta.outcome <- b[1]
          a$adj.se.outcome <- b[2]
        } else {
          message("   p->y:  ", a$SNP, "\t- ", sig$outcome[i])
          a$what <- "p->y"
          a$candidate.beta.exposure <- NA
          a$candidate.se.exposure <- NA
          a$candidate.beta.outcome <- sigo$b[sigo$id.exposure == sig$id.outcome[i]]
          a$candidate.se.outcome <- sigo$se[sigo$id.exposure == sig$id.outcome[i]]
          b <- private$bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
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
        temp$qi <- private$cochrans_q(temp$beta.outcome / temp$beta.exposure, temp$se.outcome / abs(temp$beta.exposure))
        a$adj.qi <- temp$qi[temp$SNP == a$SNP]
        a$adj.Q <- sum(temp$qi)
        
        l[[i]] <- a
      }
      
    }
    l <- bind_rows(l)
    l$d <- (l$adj.qi / l$adj.Q) / (l$qi / l$Q)
    self$output$adjustment <- l
    invisible(self$output$adjustment)
  },

############################################################################################################################
#adjustment w mvmr
  #' @description 
  #' Similar to adjusment, but when there are multiple traits associated with a single variant. 
  #' 
  #' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used.
  #' 
  #' @param tryxscan Output from \code{x$scan()}
  #' 
  #' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
  #' 
  #' @param lasso Whether to shrink the estimates of each trait within SNP. Default=TRUE.
  #' 
  #' @param proxies Look for proxies in the MVMR methods. Default = FALSE.
    adjustment.mv = function(tryxscan=self$output, lasso=TRUE, id_remove=NULL, proxies=FALSE) {
    if(!any(tryxscan$search$sig))
    {
      return(NULL)
    }
    
    l <- list()
    sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
    sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)
    sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)
    
    id.exposure <- tryxscan$dat$id.exposure[1]
    id.outcome <- tryxscan$dat$id.outcome[1]
    
    #	for each outlier SNP find its effects on 
    
    sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
    sigo1 <- subset(sig, id.outcome %in% sigo$id.exposure) %>% group_by(SNP) %>% mutate(snpcount=n()) %>% arrange(desc(snpcount), SNP)
    snplist <- unique(sigo1$SNP)
    mvo <- list()
    dat <- tryxscan$dat
    dat$qi <- private$cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
    dat$Q <- sum(dat$qi)
    dat$orig.beta.outcome <- dat$beta.outcome
    dat$orig.se.outcome <- dat$se.outcome
    dat$orig.qi <- dat$qi
    dat$outlier <- FALSE
    dat$outlier[dat$SNP %in% tryxscan$outliers] <- TRUE
    for(i in 1:length(snplist))
    {
      message("Estimating joint effects of the following trait(s) associated with ", snplist[i])
      temp <- subset(sigo1, SNP %in% snplist[i])
      candidates <- unique(temp$outcome)
      message(paste(candidates, collapse="\n"))
      if(lasso & length(candidates) == 1)
      {
        message("Only one candidate trait for SNP ", snplist[i], " so performing standard MVMR instead of LASSO")
      }
      mvexp <- suppressMessages(mv_extract_exposures(c(id.exposure, unique(temp$id.outcome)), find_proxies=proxies))
      mvout <- suppressMessages(extract_outcome_data(mvexp$SNP, id.outcome))
      mvdat <- suppressMessages(mv_harmonise_data(mvexp, mvout))
      if(lasso & length(candidates) > 1)
      {
        message("Performing shrinkage")		
        b <- glmnet::cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
        c <- coef(b, s = "lambda.min")
        keeplist <- unique(c(rownames(c)[!c[,1] == 0], tryxscan$dat$id.exposure[1]))
        message("After shrinkage keeping:")
        message(paste(keeplist, collapse="\n"))
        mvexp2 <- subset(mvexp, id.exposure %in% keeplist)
        remsnp <- group_by(mvexp2, SNP) %>% summarise(mp = min(pval.exposure)) %>% filter(mp > 5e-8) %$% as.character(SNP)
        mvexp2 <- subset(mvexp2, !SNP %in% remsnp)
        mvout2 <- subset(mvout, SNP %in% mvexp2$SNP)
        mvdat2 <- suppressMessages(mv_harmonise_data(mvexp2, mvout2))
        mvo[[i]] <- mv_multiple(mvdat2)$result
      } else {
        mvo[[i]] <- mv_multiple(mvdat)$result
      }
      mvo[[i]] <- subset(mvo[[i]], !exposure %in% tryxscan$dat$exposure[1])
      temp2 <- with(temp, tibble(SNP=SNP, exposure=outcome, snpeff=beta.outcome, snpeff.se=se.outcome, snpeff.pval=pval.outcome))
      mvo[[i]] <- merge(mvo[[i]], temp2, by="exposure")
      boo <- with(subset(dat, SNP == snplist[i]), 
                  private$bootstrap_path(
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
    
    dat$qi[dat$mr_keep] <- private$cochrans_q(dat$beta.outcome[dat$mr_keep] / dat$beta.exposure[dat$mr_keep], dat$se.outcome[dat$mr_keep] / abs(dat$beta.exposure[dat$mr_keep]))
    self$output$adjustment.mv <- list(mvo=mvo, dat=dat)
    invisible(self$output)
  },


##############################################################################################################################
#Analyse
  #' @description  
  #' This returns various heterogeneity statistics, IVW estimates for raw, adjusted and outlier removed datasets, and summary of peripheral traits detected etc.
  #'
  #' @param tryxscan Output from \code{x$scan()}.
  #' 
  #' @param plot Whether to plot or not. Default is TRUE.
  #' 
  #' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
  #' 
  #' @param duplicate_outliers_method Sometimes more than one trait will associate with a particular outlier. TRUE = only keep the trait that has the biggest influence on heterogeneity.
  analyse = function(tryxscan=self$output, plot=TRUE, id_remove=NULL, filter_duplicate_outliers=TRUE) {
    
    analysis <- list()
    x$adjustment(tryxscan=self$output, id_remove=id_remove)
    adj_full <- tryxscan$adjustment
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
    dat$qi <- private$cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
    
    analysis$Q$full_Q <- adj$Q[1]
    analysis$Q$num_reduced <- sum(adj$adj.qi < adj$qi)
    analysis$Q$mean_d <- mean(adj$d)
    
    
    temp <- subset(adj, select=c(SNP, adj.beta.exposure, adj.beta.outcome, adj.se.exposure, adj.se.outcome, candidate))
    temp$qi <- private$cochrans_q(temp$adj.beta.outcome / temp$adj.beta.exposure, temp$adj.se.outcome / abs(temp$adj.beta.exposure))
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
    dat_adj$qi <- private$cochrans_q(dat_adj$beta.outcome / dat_adj$beta.exposure, dat_adj$se.outcome / abs(dat_adj$beta.exposure))
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
    tt$qi <- private$cochrans_q(tt$beta.outcome / tt$beta.exposure, tt$se.outcome / abs(tt$beta.exposure))
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
    self$output$analysis <- analysis
    
    
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
      self$output$analysis$plot <- p
    }
    invisible(self$output)
  },


##############################################################################################################################
#Analyse w mvmr
  #' @description  
  #' Similar to analyse, but when there are multiple traits associated with a single variant. 
  #' 
  #' @param tryxscan Output from \code{x$scan()}
  #' 
  #' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
  #' 
  #' @param lasso Whether to shrink the estimates of each trait within SNP. Default=TRUE.
  #' 
  #' @param proxies Look for proxies in the MVMR methods. Default = FALSE.
  analyse.mv = function(tryxscan=self$output, lasso=TRUE, plot=TRUE, id_remove=NULL, proxies=FALSE) {
    x$adjustment.mv(tryxscan=self$output, lasso=lasso, id_remove=id_remove, proxies=proxies)
    adj <- tryxscan$adjustment.mv
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
    self$output$analyse.mv <- analysis
    invisible(self$output)
  },


############################################################################################################################
# Plots
############################################################################################################################
#' @description
#' Draw a Manhattan style plot for candidate traits-outcome/exposure associations.
#' 
#' @param what Analyse candidate-exposure ('exposure') or candidate-outcome ('outcome') associations
#' 
#' @param id_remove List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.
#' 
#' @param y_scale The scaling function to be applied to y scale.
#' 
#' @param label Display the names of the traits on the graph.
   manhattan_plot = function(what="outcome", id_remove=NULL, y_scale=NULL, label = TRUE){
      
      cpg <- require(ggplot2)
      if(!cpg)
      {
        stop("Please install the ggplot2 package")
      }
      cpg <- require(ggrepel)
      if(!cpg)
      {
        stop("Please install the ggrepel package")
      }
      
      #Open & clean data
      #mr outcome: candidate traits-outcome / candidate traits-exposure / exposure-candidate traits
      
      #X-axis preparing: 
      #Each trait needs to be shown in different phenotype groups, along the X axis 
      #Sort the traits based on their p-value for MR results
      
      stopifnot(what %in% c("exposure", "outcome"))
      dat <- self[["output"]][[paste0("candidate_", what, "_mr")]]
      temp <- subset(dat, !id.exposure %in% id_remove)
      temp <- temp[order(temp$pval_adj),]

      #Numbering row according to phenotype group where numeric variable is required.
      temp <- temp %>% mutate(id = row_number())
      
      #To compute the (mimic) cumulative position for the traits
      #pos = temp %>% 
      #      #group_by(cat) %>% 
      #      summarize(center=( max(id) + min(id) ) / 2 )
      
      #
      name <- tidyr::separate(temp, exposure, sep="\\|", c("exposure", "blank", "id"))
      temp <- subset(name, select= -c(blank))
      
      #manhattan plot
      p <- ggplot(temp, aes(x=id, y=as.numeric(-log10(pval_adj)*sign(b)))) +
        geom_point(alpha=0.8, size=3) +
        geom_point(colour = "snow", size = 1.5) +
        
        # custom X axis
        # scale_x_discrete("Phenotype", breaks=gd$center, labels = c(gd$cat))+
        #scale_x_continuous(breaks=pos$center) +
        scale_y_continuous(limits = y_scale) +
        # scale_y_continuous(expand = c(0, 0) ) +   
        
        # Add highlighted points
        geom_point(data=subset(temp, sig), size=3) +
        
        # Add line
        geom_hline(yintercept=0,colour="black", alpha=I(2/3)) +
        
        # Custom the titles
        #ggtitle("Effect of the candidate traits on the outcome") +
        # xlab("cat") + 
        ylab("-Log10(Adjusted P-value) x sign(beta)") + 
        xlab(paste("Association of the candidate traits and ", temp$outcome[1], sep=""))+
        
        # Custom the theme:
        theme_bw() +
        theme( 
          legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x=element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.title = element_text(hjust = 0.5)
          
        )     
      
      # Add label using ggrepel to avoid overlapping
      if(label)
      {
        p + geom_label_repel(data=subset(temp, sig), aes(label=exposure), size=2)
      }
      
      self$output$plots$manhattan_plot <- p
      invisible(self$output)
    }

  ))
      
   
