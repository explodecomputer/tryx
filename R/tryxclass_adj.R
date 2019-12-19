Tryx$set("public", "adjustment", function(dat= self$output$dat, tryxscan=self$output, id_remove=NULL) {
  if(!any(tryxscan$search$sig))
  {
    return(NULL)
  }
  
  l <- list()
  sig <- subset(tryxscan$search, sig & !id.outcome %in% id_remove)
  sige <- subset(tryxscan$candidate_exposure_mr, sig & !id.exposure %in% id_remove)
  sigo <- subset(tryxscan$candidate_outcome_mr, sig & !id.exposure %in% id_remove)
  
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
  self$output$adjustment <- l
  invisible(self$output$adjustment)
  }
)


Tryx$set("public", "adjustment.mv", function(dat= self$output$dat, tryxscan=self$output, lasso=TRUE, id_remove=NULL, proxies=FALSE) {
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
  dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
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
  self$output$adjustment.mv <- list(mvo=mvo, dat=dat)
  invisible(self$output)
}
)

