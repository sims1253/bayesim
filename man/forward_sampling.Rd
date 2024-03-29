% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inverse_forward_sampling.R
\name{forward_sampling}
\alias{forward_sampling}
\title{Simulate a new dataset using forward sampling.}
\usage{
forward_sampling(x, i, n, ...)
}
\arguments{
\item{x}{An object of class brmsfit or a list of brmsfit objects.}

\item{i}{The index of a single posterior draw to simulate a dataset for.
The index is passed to \code{\link{posterior_predict}}'s "draw_ids" argument.}

\item{n}{The number of samples for the newly simulated dataset.}

\item{...}{Potential additional arguments.}
}
\value{
A data.frame containing n observations for each variable in the fit.
}
\description{
Simulate a new dataset using forward sampling.
}
