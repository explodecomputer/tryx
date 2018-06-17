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

temp <- subset(info, select=c(LDHub.filename, id, trait))
sbp <- merge(sbp, temp, by.x="Trait1", by.y="LDHub.filename")
sbp <- merge(sbp, temp, by.x="Trait2", by.y="LDHub.filename")

chd <- merge(chd, temp, by.x="Trait1", by.y="LDHub.filename")
chd <- merge(chd, temp, by.x="Trait2", by.y="LDHub.filename")

id_remove <- unique(c(
	as.character(subset(sbp, abs(rg) > 0.5)$id.y),
	as.character(subset(chd, abs(rg) > 0.5)$id.y)
))

library(tryx)
ex <- extract_instruments("UKB-a:360")
ou <- extract_outcome_data(ex$SNP, 7)
d <- harmonise_data(ex, ou)
tryxscan <- tryx.scan(d)
tryxscan <- tryx.sig(tryxscan)

# Looking at the significant candidate traits, remove the following because they are too similar to SBP or CHD
id_remove <- c("UKB-a:61","UKB-a:116","UKB-a:117","UKB-a:121","UKB-a:131","UKB-a:132","UKB-a:148","UKB-a:149","UKB-a:150","UKB-a:161","UKB-a:162","UKB-a:164","UKB-a:168","UKB-a:177","UKB-a:193","UKB-a:202","UKB-a:207","UKB-a:211","UKB-a:215","UKB-a:218","UKB-a:222","UKB-a:359","UKB-a:360","UKB-a:392","UKB-a:435","UKB-a:437","UKB-a:449","UKB-a:450","UKB-a:452","UKB-a:489","UKB-a:490","UKB-a:531","472","733","UKB-a:41","UKB-a:62","UKB-a:63","UKB-a:108","UKB-a:118","UKB-a:165","UKB-a:176","UKB-a:180","UKB-a:191","UKB-a:201","UKB-a:209","UKB-a:217","UKB-a:345","UKB-a:434","UKB-a:436","UKB-a:448","UKB-a:488","6","7","8","9","113","122","229","UKB-a:532","UKB-a:533","UKB-a:534","298","UKB-a:568","760","766","1108","1109","1110","798","UKB-a:581","1092","1094","UKB-a:24","UKB-a:310","UKB-a:22")

ta <- tryx.analyse.mv(tryxscan, id_remove=id_remove)
ta$plot
ta$estimate







tryxscan <- tryx.sig(tryxscan)
ta <- tryx.analyse(tryxscan)



## Testing

n <- 10000
g <- rnorm(n)
xy <- 3
py <- 2

x <- rnorm(n) + g * 4
p <- rnorm(n) + g * 2
y <- rnorm(n) + x * xy + p * py

(gx <- lm(x ~ g)$coefficients[2])
(gp <- lm(p ~ g)$coefficients[2])
(gy <- lm(y ~ g)$coefficients[2])


(gy_x <- gy - gp * py)
(gy_p <- gy - gx * xy)




n <- 10000
px <- 1
xy <- 3
py <- 2

g <- rnorm(n)
p <- rnorm(n) + g * 2
x <- rnorm(n) + g * 4 + p * px
y <- rnorm(n) + x * xy + p * py

(gx <- lm(x ~ g)$coefficients[2])
(gp <- lm(p ~ g)$coefficients[2])
(gy <- lm(y ~ g)$coefficients[2])


(gy_x <- gy - gp * py - (gx - gp * px))
(gy_p <- gy - gx * xy)




n <- 10000
g <- rnorm(n)
xy <- 3
p1y <- 2
p2y <- 5

x <- rnorm(n) + g * 4
p1 <- rnorm(n) + g * 2
p2 <- rnorm(n) + g * 1

y <- rnorm(n) + x * xy + p1 * p1y + p2 * p2y

(gx <- lm(x ~ g)$coefficients[2])
(gp1 <- lm(p1 ~ g)$coefficients[2])
(gp2 <- lm(p2 ~ g)$coefficients[2])
(gy <- lm(y ~ g)$coefficients[2])


(gy_x <- gy - gp1 * p1y - gp2 * p2y)



n <- 10000
p1x <- 1
xy <- 3
p1y <- 2

g <- rnorm(n)
p1 <- rnorm(n) + g * 8
x <- rnorm(n) + g * 4 + p1 * p1x
y <- rnorm(n) + x * xy + p1 * p1y

(gx <- lm(x ~ g)$coefficients[2])
(gp1 <- lm(p1 ~ g)$coefficients[2])
(gy <- lm(y ~ g)$coefficients[2])


gy - gp1 * p1y

(gy_x <- gy - gp1 * p1y - (gx - gp1 * p1x))

gy_x / gx



(gy_x <- gy )


4 * 3 + 8 * 2 + 8 * 1 * 3
4 * 3 + 8 * 1 * 3

gx * xy + gp * px * xy
(gx + gp * px) * xy

(gx <- lm(x ~ g)$coefficients[2])


gy - gxh * xy - 





n <- 10000
g <- rnorm(n)
xy <- 3
p1y <- 2
p2y <- 5
p3y <- 3


p1 <- rnorm(n) + g * 2
p2 <- rnorm(n) + g * 1
p3 <- rnorm(n) + g * 6
x <- rnorm(n) + g * 4 + p3 * -2

y <- rnorm(n) + x * xy + p1 * p1y + p2 * p2y + p3 * p3y

(gx <- lm(x ~ g)$coefficients[2])
(gp1 <- lm(p1 ~ g)$coefficients[2])
(gp2 <- lm(p2 ~ g)$coefficients[2])
(gp3 <- lm(p3 ~ g)$coefficients[2])
(gy <- lm(y ~ g)$coefficients[2])


(gy_x <- gy - gp1 * p1y - gp2 * p2y - gp3 * p3y)

gy_x / gx





