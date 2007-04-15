\name{centromere}
\alias{centromere}
\title{Coordinates of centromere}
\description{
  Extracts coordinates of centromere.
}
\usage{
centromere(annotation, chromosome)
}
\arguments{
  \item{annotation}{ data.frame of chromosome annotation supplied by this package}
  \item{chromosome}{ character vector indicating which chromosomes to
    obtain centromere coordinates}
}
\details{
}
\value{
  A data.frame containing start and stop sites of the centromere for all
  chromosomes in argument chromosome.
}
\references{}
\author{Robert Scharpf}
\note{ This is a static table.  
}
\seealso{}
\examples{
data(annSnpset)
centromere(chromosomeAnnotation(annSnpset), c(1:5, "X"))
}
\keyword{manip}

