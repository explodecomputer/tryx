library(devtools)
load_all()


o <- tryx.simulate(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 2, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0, bu1x = 0.6, bu1y = 0.4, bxu3 = 0.3, bu3y = 0, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = "bonferroni", outliers_known = "all") %>% tryx.sig
tryx.analyse(o)


set.seed(1)
o1 <- tryx.simulate(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 2, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0, bu1x = 0.6, bu1y = 0.4, bxu3 = 0.3, bu3y = 0, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = "bonferroni", outliers_known = "all") %>% tryx.sig
dev.new()
tryx.analyse(o1)

set.seed(1)
o2 <- tryx.simulate(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 2, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0, bu1x = 0.6, bu1y = 0.4, bxu3 = 0.3, bu3y = 0, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = "bonferroni", outliers_known = "detected") %>% tryx.sig
tryx.analyse(o2)


set.seed(1)
o2 <- tryx.simulate(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 25, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0, bu1x = 0.6, bu1y = 0, bxu3 = 0.3, bu3y = 0, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = 0.05, outliers_known = "detected", directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o2)


set.seed(1)
o2 <- tryx.simulate(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 25, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0.2, bu1x = 0.6, bu1y = 0, bxu3 = 0.3, bu3y = 0, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = 0.05, outliers_known = "detected", directional_bias=FALSE) %>% tryx.sig
tryx.analyse(o2)



set.seed(1)
o2 <- tryx.simulate(nid = 10000, ngx = 30, ngu1 = 30, ngu2 = 30, nu2 = 20, ngu3 = 30, vgx = 0.2, vgu1 = 0.6, vgu2 = 0.2, vgu3 = 0.2, bxy = 0, bu1x = 0.6, bu1y = 0, bxu3 = 0.3, bu3y = 0.2, vgxu2 = 0.2, vu2y = 0.2, mininum_instruments = 10, instrument_threshold = "bonferroni", outlier_threshold = 0.05, outliers_known = "detected", directional_bias=FALSE) %>% tryx.sig
tryx.analyse(o2)







## old






nid, nu1, nu2, bxu3=0, bu3y=0, bxy=3,

o <- tryx.simulate(nid=10000, nu1=0, nu2=2, bxu3=0, bu3y=0, bxy=0, bu4x=0, bu4y=0) %>% tryx.sig
tryx.analyse(o)
glimpse(o)
o$search

o <- tryx.simulate(10000, 2, 2, 0.5, 0.5, 0, outlier_threshold='nominal') %>% tryx.sig
o <- tryx.simulate(10000, 2, 2, 0.5, 0.5, 0, outlier_threshold='nominal', outliers_known='detected') %>% tryx.sig
o <- tryx.simulate(10000, 2, 2, 0.5, 0.5, 0, outlier_threshold='nominal', outliers_known='all') %>% tryx.sig


o <- tryx.simulate(10000, 1, 10, 0.5, 0.5, 0, outlier_threshold='nominal', outliers_known='detected', directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o)


o <- tryx.simulate(10000, 1, 10, 0.5, 0.5, 0, outlier_threshold='nominal', outliers_known='detected', directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o)


o <- tryx.simulate(10000, 1, 10, 0.5, 0.5, 0, outlier_threshold='nominal', outliers_known='all', directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o)

o <- tryx.simulate(10000, 1, 10, 0.5, 0.5, 1, outlier_threshold='nominal', outliers_known='all', directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o)

o <- tryx.simulate(10000, 1, 10, 0.5, 0.5, 0, outlier_threshold='nominal', outliers_known='known', directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o)


o <- tryx.simulate(10000, 1, 28, 0, 0.5, 0, outlier_threshold='nominal', outliers_known='detected', directional_bias=FALSE) %>% tryx.sig
tryx.analyse(o)

o <- tryx.simulate(10000, 1, 28, 0, 0.5, 0, outlier_threshold='nominal', outliers_known='detected', directional_bias=TRUE) %>% tryx.sig
tryx.analyse(o)


tryx.adjustment(o)
glimpse(o)
o$search


o <- tryx.simulate(10000, 2, 2, 0, 0, 3)
glimpse(o)
o$search

o <- tryx.simulate(10000, 2, 2, 0, 3, 0)
glimpse(o)
o$search

o <- tryx.simulate2(10000, 2, 2, 3, 0, 0)
glimpse(o)
o$search


p <- tryx.adjustment(o)
o$search %>% sapply(., class)



o <- tryx.simulate(10000,2,2)
o <- tryx.sig(o)
p <- tryx.analyse(o)


# Simulate 100 genotypes
g <- make_geno(100000, 80, 0.5)

# Choose effect sizes for instruments for each trait
effs1 <- choose_effects(50, 0.05)
effs2 <- choose_effects(50, 0.05)

# Create X1 and X2, where they overlap some variants
x1 <- make_phen(effs1, g[,1:50])
x2 <- make_phen(effs2, g[,31:80])

# Create Y - x1 has a -0.3 influence on it and x2 has a +0.3 influence on it
y <- make_phen(c(-0.3, 0.3), cbind(x1, x2))

# Perform separate MR on each
dat1 <- get_effs(x1, y, g)
dat2 <- get_effs(x2, y, g)
library(TwoSampleMR)
mr(subset(dat1, pval.exposure < 5e-8))
mr(subset(dat2, pval.exposure < 5e-8))

# Do multivariable MR
# First get the effects for x1, x2 and y, and put them in mv format
mvdat <- make_mvdat(list(a=x1, b=x2), y, g)

# Perform MV MR
mv_multiple(mvdat)

