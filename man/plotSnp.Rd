\name{plotSnp}
\alias{plotSnp}
\alias{plotSnp,eSet-method}
\title{Plots copy number and genotype calls against physical position }
\description{

  Plots copy number and genotype call versus physical position for
  an arbitrary number of samples and chromosomes.

}
\usage{
plotSnp(object, ...)
}

\arguments{
  \item{object}{ An object extending \code{eSet}}
  \item{\dots}{ Arguments to par}
}
\details{
  See examples in the vignette
}
\value{
  
}

\references{}

\author{}

\note{
}
\seealso{\code{\link{par}}}
\examples{

data(sample.snpset)
plotSnp(sample.snpset[chromosome(sample.snpset) == "1", 1], use.chromosome.size=TRUE)
}
\keyword{hplot}

