% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcov.R
\name{vcov.oglmx.selection}
\alias{vcov.oglmx.selection}
\title{Calculate Variance-Covariance Matrix for a Fitted Model Object
Returns the variance-covariance matrix of the main parameters of a
fitted model object of class \code{oglm}.See \link[stats]{vcov} for
more details}
\usage{
\method{vcov}{oglmx.selection}(object, ...)
}
\arguments{
\item{object}{
    a fitted model object, typically.  Sometimes also a
    \code{\link{summary}()} object of such a fitted model.
  }

\item{...}{
    additional arguments for method functions.  For the
    \code{\link[stats]{glm}} method this can be used to pass a
    \code{dispersion} parameter.}
}
\description{
Calculate Variance-Covariance Matrix for a Fitted Model Object
Returns the variance-covariance matrix of the main parameters of a
fitted model object of class \code{oglm}.See \link[stats]{vcov} for
more details
}
