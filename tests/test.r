library(devtools)
load_all()
library(RadialMR)
library(dplyr)
library(parallel)
library(ggplot2)

ao <- available_outcomes()
a <- extract_instruments("UKB-a:360")
b <- extract_outcome_data(a$SNP, 7)
dat <- harmonise_data(a, b)
tryxscan <- tryx.scan(dat, mr_method="mr_ivw")
save(tryxscan, file="sbp-chd.rdata")

load("../data/sbp-chd.rdata")
ls()


ao <- available_outcomes()
info <- read.csv("../data/info.csv")



tryxscan <- tryx.sig(tryxscan)
out <- tryx.analyse(tryxscan)
out

tryxanalysis <- tryx.analyse.mv(tryxscan)

a1 <- tryx.adjustment(tryxscan)
a2 <- tryx.adjustment(tryxscan, duplicate_outliers_method = "none")
a3 <- tryx.adjustment(tryxscan, duplicate_outliers_method = "mv")
a4 <- tryx.adjustment(tryxscan, duplicate_outliers_method = "mv_lasso")



temp <- subset(outlierscan$search, pval.outcome < 5e-8)
# remove any cholesterol related traits
temp <- subset(temp, 
!grepl("simvastatin", outcome, ignore.case = TRUE) &
!grepl("atorvastatin", outcome, ignore.case = TRUE) &
!grepl("cholesterol", outcome, ignore.case = TRUE) &
!grepl("heart disease", outcome, ignore.case = TRUE) &
!grepl("HDL", outcome, ignore.case = TRUE) &
!grepl("LDL", outcome, ignore.case = TRUE) &
!grepl("myocardial", outcome, ignore.case = TRUE) &
!grepl("ischaemic", outcome, ignore.case = TRUE) &
!grepl("angina", outcome, ignore.case = TRUE)
)

temp$outcome[grepl("mass", temp$outcome, ignore.case=TRUE)]
temp$outcome[grepl("ischaemic", temp$outcome, ignore.case=TRUE)]

outlierscan <- outlier_sig(outlierscan)
adj <- outlier_adjustment(outlierscan)
plot_outliers(outlierscan, adj %>% arrange(d) %>% filter(!duplicated(SNP)))


outlier_network(outlierscan)


urate <- extract_instruments(1055)
egfr <- extract_outcome_data(urate$SNP, 1105)
chd <- extract_outcome_data(urate$SNP, 7)
urate_egfr <- harmonise_data(urate, egfr)
urate_chd <- harmonise_data(urate, chd)
outlierscan <- tryx.scan(urate_egfr, mr_method="strategy1", use_proxies=FALSE)
outlierscan <- outlier_sig(outlierscan)
volcano_plot(rbind(outlierscan$candidate_exposure_mr, outlierscan$candidate_outcome_mr))
outlier_network(outlierscan)
temp <- outlier_adjustment(outlierscan)

outlierscan2 <- tryx.scan(urate_egfr, mr_method="mr_ivw")
outlierscan2 <- outlier_sig(outlierscan2)
tem2 <- outlier_adjustment(outlierscan)

tem2$d <- (tem2$adj.qi / tem2$adj.Q) / (tem2$qi / tem2$Q)


res <- mr(dat)
ind <- dat$beta.exposure < 0
dat$beta.exposure[ind] <- abs(dat$beta.exposure[ind])
dat$beta.outcome[ind] <- dat$beta.outcome[ind] * -1
mr_scatter_plot(res, dat)[[1]] +
geom_point(data=subset(dat, SNP %in% l$SNP), colour="red") +
geom_point()


plot_outliers(outlierscan, adj %>% arrange(d) %>% filter(!duplicated(SNP)))

# is the outlier adjustment doing the right thing?
dev.new()
RadialMR(outlierscan$dat$beta.exposure,outlierscan$dat$beta.outcome,outlierscan$dat$se.exposure,outlierscan$dat$se.outcome,outlierscan$dat$SNP,"BOTH","YES","NO",0.05,"NO")$plot


a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 7)
dat <- harmonise_data(a, b)
outlierscan <- outlier_scan(dat, mr_method="mr_ivw")
outlierscan <- outlier_sig(outlierscan)
adj <- outlier_adjustment(outlierscan)
plot_outliers(outlierscan, adj %>% arrange(d) %>% filter(!duplicated(SNP)))

