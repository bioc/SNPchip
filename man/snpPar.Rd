\name{snpPar}
\alias{snpPar}
\alias{snpPar<-}
\title{Accessor for graphical parameters}
\description{
Accessor for the list of graphical parameters in objects of class \code{ParESet}.
}
\usage{
snpPar(object)
}
\arguments{
  \item{object}{object of class \code{ParESet}}
}
\value{
  \code{list}
}
\author{RS}
\seealso{ \code{\link{par}}, \code{\link[graphics]{layout}}}
\examples{
data(sample.snpset)
object <- plotSnp(sample.snpset)
str(snpPar(object))
\dontrun{
show(object)
}
}
\keyword{methods}

