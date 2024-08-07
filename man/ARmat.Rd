% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_test.R
\name{ARmat}
\alias{ARmat}
\title{Create Covariance Matrix with Autoregressive Correlation Structure}
\usage{
ARmat(time, sigma, rho)
}
\arguments{
\item{time}{vector of time-points.}

\item{sigma}{vector of the standard deviations assumed at the different time point.}

\item{rho}{correlation between adjacent time-points of distance 1, see details.}
}
\value{
A \eqn{k} by \eqn{k} covariance matrix of AR(1) structure, with \eqn{k} being the number of the time-points.
}
\description{
Creates a covariance matrix of structure AR(1) (autoregressive of order 1) based on given time points, correlation between adjacent time points and standard deviations of each time point.
}
\details{
The correlation between two time points \eqn{a} and \eqn{b} and AR(1) correlation parameter \eqn{\rho} is \eqn{\rho^{|a-b|}}.
}
\examples{
ARmat(1:3,rep(1,3),0.5)
K<-ARmat(1:3,c(0.5,0.7,1),0.5)
cov2cor(K)

}
\author{
Robin Ristl \email{robin.ristl@meduniwien.ac.at}
}
