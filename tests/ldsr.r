library(devtools)
load_all()

library(RadialMR)
library(dplyr)
library(parallel)
library(ggplot2)


load("data/sbp-chd.rdata")

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

id_remove <- c(id_remove, "UKB-a:24", "UKB-a:310", "UKB-a:22")

ta <- tryx.analyse(tryxscan, id_remove=id_remove)
ta$plot









tryxscan <- tryx.sig(tryxscan)
ta <- tryx.analyse(tryxscan)



## Testing

n <- 10000
g <- rnorm(n)
x <- rnorm(n) + g * 4
p <- rnorm(n) + g * 2
y <- rnorm(n) + x * 3 + p * 2

gx <- lm(x ~ g)$coefficients[2]
gp <- lm(x ~ p)$coefficients[2]
gy <- lm(y ~ g)$coefficients[2]

xy <- 3
py <- 2

gy_x <- gy - gp * py
gy_p <- gy - gx * xy




n <- 10000
g <- rnorm(n)
p <- rnorm(n) + g * 2
x <- rnorm(n) + g * 4 + p
y <- rnorm(n) + x * 3 + p * 2

gx <- lm(x ~ g)$coefficients[2]
gp <- lm(x ~ p)$coefficients[2]
gy <- lm(y ~ g)$coefficients[2]

xy <- 3
py <- 2
px <- 1

(gy_x <- gy - gp * py - (gx - gp * px) * xy)


(gy_p <- gy - gx * xy)



