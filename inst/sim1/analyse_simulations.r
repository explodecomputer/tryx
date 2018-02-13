library(tidyverse)
load("sim3.rdata")

temp <- lapply(1:length(l), function(x){
	a <- l[[x]]$estimates
	a$sim <- x
	return(a)
}) %>% bind_rows %>% inner_join(param)

temp$nu <- temp$nu1 + temp$nu2

dim(temp)
table(temp$est)
temp$diff <- temp$b - temp$bxy
ggplot(temp, aes(y=diff, x=est)) +
geom_boxplot() +
facet_grid(nid ~ bxy)

temp2 <- group_by(temp, est, bxy, nid, nu) %>%
summarise(pow=sum(pval < 1e-5)/n(), isq=mean(Isq))

ggplot(temp2, aes(y=pow, x=nu)) +
geom_line(aes(colour=est)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual")

temp3 <- group_by(temp, est, bxy, nid) %>%
summarise(pow=sum(pval < 1e-5)/n(), isq=mean(Isq))

ggplot(temp3, aes(y=pow, x=est)) +
geom_point() +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90))

ggplot(temp3, aes(y=isq, x=est)) +
geom_point() +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90))