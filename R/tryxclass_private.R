
Tryx$set("private", "cochrans_q", function(b, se) {
  xw <- sum(b / se^2) / sum(1/se^2)
  qi <- (1/se^2) * (b - xw)^2
  return(qi)
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
