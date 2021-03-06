% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additional_functions.R
\name{Laplace_approx}
\alias{Laplace_approx}
\title{Laplace_approx}
\usage{
Laplace_approx(u, parameters, h = 1e-04)
}
\arguments{
\item{u}{A number}

\item{parameters}{A list with given parameters to function:
\code{k1m}, \code{k2m}, \code{A}, \code{cb}, \code{tau}, \code{func_mu},
\code{func_sigma}, \code{x}, \code{theta0}}

\item{h}{numerical derivative parameter}
}
\value{
value of Laplace_approx function
}
\description{
\eqn{\mathcal{L}(\overline{\mathbf{q}_n},\overline{u})} for the \eqn{\tau} numerator
}
