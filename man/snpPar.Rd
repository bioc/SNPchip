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
\details{
}
\value{
  \code{list}
}
\references{ }
\author{}
\note{ 
}
\seealso{ \code{\link{par}}, \code{\link[graphics]{layout}}}
\examples{
data(sample.snpset)
gp <- new("ParSnpSet")
object <- getPar(gp, sample.snpset)
str(snpPar(object))

##or equivalently
show(object)
}
\keyword{methods}

