\name{plotSnp}
\alias{plotSnp}
\alias{plotSnp,ParSnpCopyNumberSet,SnpCopyNumberSet-method}
\alias{plotSnp,ParSnpCallSet,SnpCallSet-method}
\alias{plotSnp,ParSnpSet,oligoSnpSet-method}
\alias{plotSnp,SnpLevelSet-method}
\title{Plots copy number and genotype calls against physical position }
\description{

  Returns an object inheriting from class \code{ParESet} -- essentially,
  a list of default graphical parameters that can be modified as needed.
  The \code{show} method for this class plots the copy number and
  genotype calls versus physical position for an arbitrary number of
  samples and chromosomes.

}
\usage{
plotSnp(object, hmmPredict, ...)
}

\arguments{
  \item{object}{ An object extending \code{ParESet}}
  \item{hmmPredict}{An object of class \code{HmmPredict}}
  \item{\dots}{ Additional arguments to the initialization methods of
  the \code{ParESet} classes}
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
graph.par <- plotSnp(chr1)
class(graph.par)
show(graph.par) ##makes the plot
}
\keyword{hplot}

