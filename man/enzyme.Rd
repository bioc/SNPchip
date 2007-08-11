\name{enzyme}
\alias{enzyme}
\alias{enzyme,eSet-method}
\alias{enzyme,AnnotatedDataFrame-method}
\title{Enzyme used to digest DNA}
\description{When two chips have been combined (e.g., Xba and Hind),
  enzyme is an accessor for whether a particular SNP was interrogated on
  the Xba or Hind chip.  A better name for this accessor might be
  'chip', suggestions welcome.
}
\usage{
enzyme(object)
}
\arguments{
  \item{object}{one of the classes defined in \code{SNPchip}}
}
\details{
}
\value{
  Returns the enzyme for each indicated SNP.
}
\author{Robert Scharpf}
\note{
}
\seealso{}
\examples{
}
\keyword{methods}
