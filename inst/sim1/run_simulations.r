library(tryx)
library(parallel)
library(tidyverse)

param <- expand.grid(
	nid = c(5000, 10000),
	bxy = c(-0.5, -0.2, 0.2, 0.5),
	nu1 = c(1, 5, 10, 14),
	nu2 = c(1, 5, 10, 14),
	outliers_known = c(TRUE),
	simr = c(1:100)
)
param$sim <- 1:nrow(param)



l <- mclapply(1:nrow(param), function(i)
{
	set.seed(i)
	out <- simulate.tryx(param$nid[i], param$nu1[i], param$nu2[i], param$bxy[i], outliers_known=param$outliers_known[i])
	out <- outlier_sig(out)
	return(tryx.analyse(out, plot=FALSE))
}, mc.cores=16)

save(l, param, file="sim.rdata")
