% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigma.R
\name{sigma.oglmx}
\alias{sigma.oglmx}
\title{Recover residual standard deviation from ordered discrete model}
\usage{
\method{sigma}{oglmx}(object, ...)
}
\arguments{
\item{object}{A \code{oglmx} model}

\item{...}{Additional arguments. Consider in particular adding
\code{newdata} to parameters}
}
\value{
Residual estimated standard deviation in vector form. With an
homoskedastic model, all values are equal
}
\description{
Extract the estimated standard deviation of the errors,
the “residual standard deviation” (misnomed also “residual standard error”).
}
\details{
This function transforms the linear model for
the standard deviation into the standard deviation
}
