% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{espls}
\alias{espls}
\title{Kernel M-Estimator for Location Scale model}
\usage{
espls(Yt, St, s, initial.values, bandwidth = 1.06 * sqrt(var(St)) *
  length(St)^(-1/5), int.of.par = c(0, 1), print = F)
}
\arguments{
\item{Yt}{parmeter of a function which is not to be optimized, usually \eqn{Y_t}}

\item{St}{regresor parameter can be X's or lag(Y_t)}

\item{s}{points at which function should be estimated}

\item{initial.values}{initial value of optimisible parameter might be a vector}

\item{bandwidth}{bandwith should be used}

\item{int.of.par}{initial parameters}

\item{print}{print during fitting}
}
\value{
Estimated location scale function at \code{s} points
}
\description{
Estimates parameters for location scale model using Kernel M-Estimator using R
optim function
}
