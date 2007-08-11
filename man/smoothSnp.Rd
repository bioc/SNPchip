\name{smoothSnp}
\alias{smoothSnp}
\alias{smoothSnp,oligoSnpSet-method}
\title{A simple nonparametric smoother for genotype and copy number}
\description{
  A loess smoother for finding regions of reduced copy number and loss
  of heterozygosity.
}
\usage{
smoothSnp(object, ...)
}
\arguments{
  \item{object}{Object derived from \code{AnnotatedSnpSet}}
  \item{\dots}{For additional arguments see details.}
}
\details{
  Additional arguments include the chromosomes (character vector) and
  samples (numerical vector) to smooth.  Options for smoothing are set
  by specifying the span and the method.  The loess method is the
  default and the only method currently available.
}

\value{ An object of the same class, e.g., \code{AnnotatedSnpSet}, but
  assayData elements copyNumber and calls are replaced by the smoothed
  values.  }

\references{ }
\author{Robert Scharpf}
\note{ 
}
\seealso{}
\examples{
}
\keyword{methods}

