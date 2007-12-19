\name{centromere}
\alias{centromere}
\title{Coordinates of centromere}
\description{
  Extracts coordinates of centromere for a particular chromosome
}
\usage{
centromere(chromosome)
}
\arguments{
  \item{chromosome}{character string: "1", ..., "22", "X", or "Y"}
}
\value{
  integer: stat and stop coordinates of centromere in basepairs
}
\author{R. Scharpf}
\examples{
centromere("1")
}
\keyword{manip}

