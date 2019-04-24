library(devtools)
load_all()

o <- tryx.simulate(10000, 2, 2, 3, 3, 0) %>% tryx.sig
tryx.adjustment(o)
tryx.analyse(o)
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

