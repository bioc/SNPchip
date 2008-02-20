\name{getPar}
\alias{getPar}
\title{Adds graphical parameters for plotting SNP data}
\description{
  Adds graphical parameters for plotting SNP data to one of the
  graphical parameter classes (e.g., \code{ParSnpSet},
  \code{ParSnpCopyNumberSet}, etc)
}

\usage{
  getPar(object, snpset, ...)
}

\arguments{
  \item{object}{ an object inheriting from \code{ParESet}}
  \item{snpset}{an object inheriting from \code{SnpLevelSet}}
  \item{\dots}{additional arguments to \code{par}}
}
\details{

Adds graphical parameters to an object inherited from class
      \code{ParESet} that depend on the data class object (an object
      inherited from class \code{SnpLevelSet}). For instance, graphical
      parameters specifying the layout depend on the number of samples
      and chromosomes in the object \code{snpset}
      
}

\value{
  An object of the same class as \code{object}
}

\author{R. Scharpf}
\seealso{\code{\link{ParESet-class}}, \code{\link{ParSnpSet-class}},
\code{\link{ParSnpCopyNumberSet-class}}, \code{\link{ParSnpCallSet-class}}}

\examples{
data(sample.snpset)
gp <- new("ParSnpSet")
##Default layout parameters:
gp$mat
##update layout parameters, etc, using information about the data object
gp <- getPar(gp, sample.snpset)
gp$mat
}

\keyword{manip}
\keyword{methods}
