\name{plotSnp}
\alias{plotSnp}
\title{Plot genome-wide or chromosome-specific views of copy number
  and/or genotype calls versus physical position.}
\description{
  Function for genome-wide plots of copy number and genotype calls.
}
\usage{
plotSnp(object, chromosomes, samples, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ \code{AnnotatedSnpCopyNumberSet} or
    \code{AnnotatedSnpSet}}
  \item{chromosomes}{vector of integers indicating which chromosomes to
    plot}
  \item{samples}{vector of integers indicating which samples to plot (or
    which columns in the assayData matrices.}
  \item{\dots}{additional arguments}
}
\details{
}
\value{
}
\author{Robert Scharpf}
\note{ 
}
\seealso{ \code{\link{plotChromosome}} }
\examples{
data(annSnpset)
##for formatting the plot layout, see the vignette
plotSnp(chromosomes = 1:10, object = annSnpset, samples = 1:2)
}
\keyword{methods}
