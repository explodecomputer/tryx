% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustment.r
\name{tryx.adjustment}
\alias{tryx.adjustment}
\title{Outlier adjustment estimation}
\usage{
tryx.adjustment(tryxscan, id_remove = NULL)
}
\arguments{
\item{tryxscan}{Output from \code{tryx.scan}}

\item{id_remove}{List of IDs to exclude from the adjustment analysis. It is possible that in the outlier search a candidate trait will come up which is essentially just a surrogate for the outcome trait (e.g. if you are analysing coronary heart disease as the outcome then a variable related to heart disease medication might come up as a candidate trait). Adjusting for a trait which is essentially the same as the outcome will erroneously nullify the result, so visually inspect the candidate trait list and remove those that are inappropriate.}
}
\value{
data frame of adjusted effect estimates and heterogeneity stats
}
\description{
How much of the heterogeneity due to the outlier can be explained by alternative pathways?
}
