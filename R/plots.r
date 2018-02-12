#' Plot volcano plot of many MR analyses
#' 
#' @param res Dataframe similar to output from mr() function, requiring id.exposure, id.outcome. Ideally only provide one MR estimate for each exposure-outcome hypothesis. Also provide a sig column of TRUE/FALSE to determine if that association is to be labelled
#' @param what Whether to plot many exposures against few outcomes (default: exposure) or few exposures against many outcomes (outcome). If e.g. 'exposure' and there are multiple outcomes then will facet by outcome
#' @export
#' @return ggplot of volcano plots
volcano_plot <- function(res, what="exposure")
{
	cpg <- require(ggplot2)
	if(!cpg)
	{
		stop("Please install the ggplot2 package")
	}
	cpg <- require(ggrepel)
	if(!cpg)
	{
		stop("Please install the ggrepel package")
	}

	stopifnot(all(c("outcome", "exposure", "b", "se", "pval") %in% names(res)))
	if(!"sig" %in% names(res))
	{
		warning("Significant associations not defined. Defaulting to p < 0.05 cutoff")
		res$sig <- res$pval < 0.05
	}

	if(what == "exposure")
	{
		form <- as.formula(". ~ outcome")
		nout <- length(unique(res$outcome))
		if(nout > 5)
			warning(nout, " might be too many outcomes to plot clearly")

	} else if(what == "outcome"){
		form <- as.formula(". ~ exposure")
		nout <- length(unique(res$exposure))
		if(nout > 5)
			warning(nout, " might be too many outcomes to plot clearly")
	} else {
		stop("'what' argument should be 'exposure' or 'outcome'")
	}

	ggplot(res, aes(x=b, y=-log10(pval))) +
	geom_vline(xintercept=0, linetype="dotted") +
	geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se)) +
	geom_point(aes(colour=sig)) +
	facet_grid(form, scale="free") +
	geom_label_repel(data=subset(res, sig), aes(label=exposure), colour="black", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
	# geom_text_repel(data=subset(res, sig), aes(label=outcome, colour=category)) +
	geom_point(aes(colour=sig)) +
	theme_bw() + 
	theme(panel.grid.major = element_blank(), 
	      panel.grid.minor = element_blank(),
	      strip.text.x=element_text(size=16)
	) +
	# scale_colour_brewer(type="qual", palette="Dark2") +
	# scale_fill_brewer(type="qual", palette="Dark2") +
	labs(x="MR effect size")
}


#' Plot results from outlier_scan in a network
#' 
#' Creates a simple network depicting the connections between outlier instruments, the original exposure and outcome traits, and the detected candidate associations
#' 
#' @param tryxscan Output from outlier_scan function
#' @export
#' @return Prints plot, and returns dataframe of the connections
tryx.network <- function(tryxscan)
{

	a <- require(igraph)
	if(!a)
	{
		stop("Please install the igraph R package")
	}
	a <- require(dplyr)
	if(!a)
	{
		stop("Please install the dplyr R package")
	}

	stopifnot("candidate_outcome_mr" %in% names(tryxscan))
	stopifnot("candidate_exposure_mr" %in% names(tryxscan))

	if(!"sig" %in% names(tryxscan$candidate_outcome_mr) | !"sig" %in% names(tryxscan$candidate_exposure_mr))
	{
		stop("Significant hits have not been specified. Please run tryx.sig.")
	}


	ao <- available_outcomes()
	grp <- rbind(
		# Instruments for exposure
		data.frame(
			from=tryxscan$outliers,
			to=ao$trait[ao$id == tryxscan$dat$id.exposure[1]],
			what="Outlier instruments"
		),
		# Outlier instruments associated with candidate traits
		data.frame(
			from=subset(tryxscan$search, sig)$SNP,
			to=subset(tryxscan$search, sig)$originalname.outcome,
			what="Search associations"
		),
		# candidate traits associated with exposure
		data.frame(
			from=ao$trait[ao$id %in% subset(tryxscan$candidate_exposure_mr, sig)$id.exposure],
			to=ao$trait[ao$id == tryxscan$dat$id.exposure[1]],
			what="Candidate traits to exposure"
		),
		data.frame(
			from=ao$trait[ao$id == tryxscan$dat$id.exposure[1]],
			to=ao$trait[ao$id %in% subset(tryxscan$exposure_candidate_mr, sig)$id.exposure],
			what="Exposure to candidate traits"
		),
		# candidate traits associated with outcome
		data.frame(
			from=ao$trait[ao$id %in% subset(tryxscan$candidate_outcome_mr, sig)$id.exposure],
			to=ao$trait[ao$id == tryxscan$dat$id.outcome[1]],
			what="Candidate traits to outcome"
		),
		data.frame(
			from=ao$trait[ao$id == tryxscan$dat$id.exposure[1]],
			to=ao$trait[ao$id == tryxscan$dat$id.outcome[1]],
			what="Main hypothesis"
		)
	)


	l <- list()
	i <- 1
	temp <- ao$id %in% subset(tryxscan$candidate_exposure_mr, sig)$id.exposure & ao$id %in% subset(tryxscan$candidate_outcome_mr, sig)$id.exposure
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="both"
		)
		i <- i + 1
	}
	temp <- ao$id %in% subset(tryxscan$candidate_exposure_mr, sig)$id.exposure
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="exposure"
		)
		i <- i + 1
	}
	temp <- ao$id %in% subset(tryxscan$candidate_outcome_mr, sig)$id.exposure
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="outcome"
		)
		i <- i + 1
	}
	temp <- ao$id %in% subset(tryxscan$search, sig)$id.outcome
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="search"
		)
		i <- i + 1
	}
	sig <- bind_rows(l) %>% filter(!duplicated(id))


	nodes <- rbind(
		data_frame(
			name=c(
				ao$trait[ao$id %in% tryxscan$dat$id.exposure[1]],
				ao$trait[ao$id %in% tryxscan$dat$id.outcome[1]]),
			id=c(tryxscan$dat$id.exposure[1], tryxscan$dat$id.outcome[1]),
			what=c("original")
		),
		data_frame(
			name=unique(tryxscan$outliers),
			id=NA,
			what="Outlier instruments"
		),
		sig
	)

	nodes$size <- 15
	nodes$size[nodes$what != "original"] <- 3
	nodes$label <- nodes$name
	nodes$label[! nodes$what %in% c("original")] <- NA
	nodes$colour <- "#386cb0"
	nodes$colour[nodes$what == "original"] <- "white"
	nodes$colour[nodes$what == "Outlier instruments"] <- "#f0027f"


	# Make layout
	layoutd <- group_by(nodes, what) %>%
		do({
			x <- .
			a <- as.data.frame(t(combn(as.character(x$name), 2)))
			# a <- data.frame(from=x$name[1], to=x$name[-1])
			names(a) <- c("from", "to")
			a$what <- x$what[1]
			a
		})

	layoutd2 <- rbind(as.data.frame(layoutd), 
		data.frame(
			from=subset(nodes, what != "original")$name,
			to=sample(subset(nodes, what == "original")$name, length(subset(nodes, what != "original")$name), replace=TRUE),
			what="temp"
		)
	)

	# Different colours for different edge types
	# Original should be large
	# original nodes should be large
	# no labels for confounders

	grp$colour <- "black"
	grp$colour[grp$what == "Outlier instruments"] <- "#7fc97f"
	grp$colour[grp$what == "Search associations"] <- "#beaed4"
	grp$colour[grp$what == "Candidate traits to outcome"] <- "#fdc086"
	grp$colour[grp$what == "Candidate traits to exposure"] <- "#bf5b17"
	grp$size <- 0.1
	grp$size[grp$what == "Main hypothesis"] <- 0.5

	layoutg <- graph_from_data_frame(layoutd2, vertices=nodes)
	l <- layout_with_fr(layoutg)
	grl <- graph_from_data_frame(grp, directed=TRUE, vertices=nodes)
	plot(grl, 
		layout=l, 
		vertex.size=V(grl)$size, 
		vertex.label=V(grl)$label,
		# edge.arrow.size=E(grl)$size, 
		edge.arrow.size=0.3,
		edge.color=E(grl)$colour,
		vertex.label.cex=0.5, 
		vertex.label.family="sans", 
		vertex.label.color="black", 
		vertex.color=V(grl)$colour,
		edge.color="red"
	)
}

