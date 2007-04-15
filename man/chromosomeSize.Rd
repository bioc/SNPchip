\name{chromosomeSize}
\alias{chromosomeSize}
\title{Get size of chromosome in number of base pairs.}
\description{
  Size of chromosome in base pairs.
}
\usage{
chromosomeSize(annotation, chromosome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation}{ data.frame contiaining chromosome-level
    annotation.}
  \item{chromosome}{Vector of character strings.}
}
\details{
}
\value{
  Numerical vector.
}
\references{}
\author{Robert Scharpf}
\note{
}
\seealso{}
\examples{
data(chromosomeAnnotation)
chromosomeSize(chromosomeAnnotation, c(1:22, "X", "Y"))
}
\keyword{manip}

