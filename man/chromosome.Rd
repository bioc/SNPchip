\name{chromosome}
\alias{chromosome}
\title{accessor for chromosome SNP annotation}
\description{
}
\usage{
chromosome(object)
}
\arguments{
  \item{object}{a \code{SnpSet} object}
}
\details{
}
\value{
  The chromosome for each SNP (a vector of strings)
}
\author{Robert Scharpf}
\note{
}
\seealso{}
\examples{
data(annSnpset)
##Get the chromosome information for the first 5 SNPs.
chromosome(annSnpset)[1:5]
}
\keyword{methods}
