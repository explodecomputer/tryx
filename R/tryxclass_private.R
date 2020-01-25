
Tryx$set("private", "strategy1", function(dat, het_threshold=0.05, ivw_max_snp=1) {
  message("First pass: running ", length(unique(paste(dat$id.exposure, dat$id.outcome))), " analyses with Wald ratio or IVW")
  a <- suppressMessages(mr(dat, method_list=c("mr_ivw", "mr_wald_ratio")))
  b <- subset(a, method == "Wald ratio")
  message(nrow(b), " analyses can only run Wald ratio")
  c <- subset(a, method == "Inverse variance weighted")
  message(nrow(c), " analyses can run IVW or mode")
  if(nrow(c) > 0)
  {
    d <- suppressMessages(mr_heterogeneity(dat, method_list="mr_ivw"))
    rerun <- subset(d, Q_pval < het_threshold & Q_df >= (ivw_max_snp-1))
    rerun <- paste(rerun$id.exposure, rerun$id.outcome)
    if(length(rerun) < 0)
    {
      message("All eligible IVW results have low heterogeneity")
      return(list(res=a, all=NULL, heterogeneity=d))
    }
    e <- subset(c, !paste(id.exposure, id.outcome) %in% rerun)
    message("Rerunning MR with weighted modal estimator for ", length(rerun), " analysis - this may be quite slow")
    f <- suppressMessages(mr(
      subset(dat, paste(id.exposure, id.outcome) %in% rerun),
      method_list="mr_weighted_mode"
    ))
    res <- rbind(b, e, f)
    return(list(res=res, all=a, heterogeneity=d))
  } else {
    return(list(res=a, all=NULL, heterogeneity=NULL))
  }
  self$output$strategy1 <- res
  invisible(self$output)
}
)


Tryx$set("private", "cochrans_q", function(b, se) {
  xw <- sum(b / se^2) / sum(1/se^2)
  qi <- (1/se^2) * (b - xw)^2
  return(qi)
  invisible(x$ouput)
}
)


Tryx$set("private", "bootstrap_path1", function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000) {
  res <- rnorm(nboot, gx, gx.se) - rnorm(nboot, gp, gp.se) * rnorm(nboot, px, px.se)
  pe <- gx - gp * px
  return(c(pe, sd(res)))
}
)


Tryx$set("private", "bootstrap_path", function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000) {
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
)


Tryx$set("private", "radialmr", function(dat, outlier=NULL) {
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
)
