library(devtools)
load_all()

library(RadialMR)
library(dplyr)
library(parallel)
library(ggplot2)


load("R/sbp-chd.rdata")

info <- read.csv("data/info.csv")
sbp <- read.csv("data/bp.csv")
chd <- read.csv("data/chd.csv")

table(chd$Trait1 %in% info$LDHub.filename)
table(sbp$Trait1 %in% info$LDHub.filename)

ao <- available_outcomes()

temp <- subset(info, select=c(LDHub.filename, id, trait))
sbp <- merge(sbp, temp, by.x="Trait1", by.y="LDHub.filename")
sbp <- merge(sbp, temp, by.x="Trait2", by.y="LDHub.filename")

chd <- merge(chd, temp, by.x="Trait1", by.y="LDHub.filename")
chd <- merge(chd, temp, by.x="Trait2", by.y="LDHub.filename")

tryxscan <- tryx.sig(tryxscan)

id_remove <- unique(c(
	as.character(subset(sbp, abs(rg) > 0.5)$id.y),
	as.character(subset(chd, abs(rg) > 0.5)$id.y)
))

id_remove <- c(id_remove, "id:UKB-a:24", "id:UKB-a:310", "id:UKB-a:24")

ta <- tryx.analyse(tryxscan, id_remove=id_remove)
ta$plot

