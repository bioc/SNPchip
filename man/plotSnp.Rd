\name{plotSnp}
\alias{plotSnp}
\alias{plotSnp,ParSnpCopyNumberSet,SnpCopyNumberSet-method}
\alias{plotSnp,ParSnpCallSet,SnpCallSet-method}
\alias{plotSnp,ParSnpSet,oligoSnpSet-method}
\title{Plots copy number and genotype calls against physical position }
\description{

  Plots copy number and genotype call versus physical position for
  an arbitrary number of samples and chromosomes.

}
\usage{
plotSnp(object, snpset, ...)
}

\arguments{
  \item{object}{ An object extending \code{ParESet}}
  \item{snpset}{An object extending \code{SnpLevelSet}}
  \item{\dots}{ Arguments to par}
}
\details{
  See examples in the vignette
}

\author{R. Scharpf}

\seealso{
  \code{\link{par}},
  \code{\link{ParESet-class}},
  \code{\link{ParSnpCallSet-class}},
  \code{\link{ParSnpCopyNumberSet-class}},
  \code{\link{ParSnpSet-class}}
}

\examples{

data(sample.snpset)
chr1 <- sample.snpset[chromosome(sample.snpset) == "1", 1]
graph.par <- new("ParSnpSet")
graph.par <- getPar(graph.par, chr1, use.chromosome.size=TRUE)
plotSnp(graph.par, chr1)
}
\keyword{hplot}

