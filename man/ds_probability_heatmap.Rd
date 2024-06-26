% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ds_probability_heatmap.R
\name{ds_probability_heatmap}
\alias{ds_probability_heatmap}
\title{Heatmaps of class probabilities}
\usage{
ds_probability_heatmap(out)
}
\arguments{
\item{out}{output of the function \code{ms_identity_map}. # TODOELI: here as well, some param name more informative of what needs to be provided?}
}
\value{
A ggplot object, with the heatmap displaying the probability of each
cell to belong to the identity classes according the model.
}
\description{
This function identifies true label groups between reference groups and clusters.
}
\examples{
#not run
heatmap<-ds_probability_heatmap(out)
}
