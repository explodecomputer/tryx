library(tidyverse)
load("result.rdata")

# Combine the results in l with the simulation parameters in param
temp <- lapply(1:length(l), function(x){
	a <- l[[x]]$estimates
	a$sim <- x
	a$no_outlier_flag <- l[[x]]$detection$no_outlier_flag
	return(a)
}) %>% bind_rows %>% inner_join(param) %>% filter(!is.na(est))
temp$no_outlier_flag[is.na(temp$no_outlier_flag)] <- FALSE
temp <- group_by(temp, sim) %>%
	do({
		x <- .
		if(x$no_outlier_flag[1])
		{
			x <- rbind(
				subset(x, est == "Raw"),
				subset(x, est == "Raw"),
				subset(x, est == "Raw")
			)
			x$est <- c("Outliers removed", "Raw", "Outliers adjusted")
		}
		x
	})

temp$nu <- temp$nu1 + temp$nu2

dim(temp)
table(temp$est)

# Rename the method
temp$outliers_known2 <- "Pleiotropy known"
temp$outliers_known2[!temp$outliers_known] <- "Pleiotropy detected"
temp$method <- paste0(temp$est, " (", temp$outliers_known2, ")")
temp$method[temp$est == "Raw"] <- "Raw"

# how different from the simulated causal effect is the estimated causal effect?
temp$diff <- temp$b - temp$bxy
temp$pdiff <- pt(abs(temp$diff)/temp$se, temp$nsnp, lower.tail=FALSE)

## Main plots

temp2 <- filter(temp, bxy >= 0) %>%
	group_by(method, nu, bxy, nid) %>%
	summarise(
		n=n(),
		psig=sum(pval < 0.05) / n(),
		bias=sum(pdiff < 0.05) / n(),
		isq=mean(Isq)
	)

ggplot(temp2, aes(y=psig, x=nu)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method")
ggsave("images/sig1.pdf")

ggplot(temp2, aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method")
ggsave("images/bias1.pdf")


temp3 <- filter(temp, bxy >= 0) %>%
	group_by(method, nu, bxy) %>%
	summarise(pow=sum(pval < 0.05)/n())
ggplot(temp3, aes(y=pow, x=nu)) +
geom_point(aes(colour=method)) +
geom_line(aes(colour=method)) +
facet_grid(. ~ bxy) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="")
ggsave("images/sig2.pdf")


temp3 <- filter(temp, bxy >= 0) %>%
	group_by(method, nu) %>%
	summarise(bias=sum(pdiff < 0.05)/n())
ggplot(temp3, aes(y=bias, x=nu)) +
geom_point(aes(colour=method)) +
geom_line(aes(colour=method)) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method")
ggsave("images/bias2.pdf")


temp3 <- filter(temp, bxy >= 0) %>%
	group_by(method, bxy) %>%
	summarise(pow=sum(pval < 0.05)/n())
ggplot(temp3, aes(y=pow, x=method)) +
geom_bar(aes(fill=as.factor(bxy)), stat="identity", position="dodge") +
scale_fill_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="", y="Proportion of estimates with p < 0.05", fill="Simulated\ncausal effect")
ggsave("images/sig3.pdf")


temp3 <- filter(temp, bxy >= 0) %>%
	group_by(method) %>%
	summarise(bias=sum(pdiff < 0.05)/n())
ggplot(temp3, aes(y=bias, x=method)) +
geom_bar(aes(fill=method), stat="identity") +
scale_fill_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="none") 
ggsave("images/bias3.pdf")





## Other plots

ggplot(temp, aes(y=diff, x=method)) +
geom_boxplot() +
facet_grid(nid ~ bxy) +
theme(axis.text.x=element_text(angle=90))

ggplot(temp, aes(y=pdiff, x=method)) +
geom_boxplot() +
facet_grid(nid ~ bxy) +
theme(axis.text.x=element_text(angle=90))

ggplot(temp, aes(y=-log10(pdiff), x=method)) +
geom_boxplot() +
facet_grid(nid ~ bxy) +
theme(axis.text.x=element_text(angle=90))



temp2 <- group_by(temp, method, bxy, nid, nu) %>%
summarise(pow=sum(pval < 1e-5)/n(), isq=mean(Isq))

ggplot(temp2, aes(y=pow, x=nu)) +
geom_line(aes(colour=method)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual")

temp3 <- group_by(temp, method, bxy, nid) %>%
summarise(pow=sum(pval < 0.05)/n(), isq=mean(Isq))

ggplot(temp3, aes(y=pow, x=method)) +
geom_point() +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
geom_hline(yintercept = 0.05, linetype="dotted") +
theme(axis.text.x=element_text(angle=90))

ggplot(temp3, aes(y=isq, x=method)) +
geom_point() +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90))


