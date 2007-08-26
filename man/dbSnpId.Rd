\name{dbSnpId}
\alias{dbSnpId}
\alias{dbSnpId,AnnotatedDataFrame-method}
\alias{dbSnpId,SnpLevelSet-method}
\title{The dbSNP identifier }
\description{
  Unique identifier for each SNP assigned by NIH
}
\usage{
dbSnpId(object)
}
\arguments{
  \item{object}{Object possibly derived from \code{eSet} or
    \code{AnnotatedDataFrame} class.}
}
\details{
}
\value{
  A (typically very long) character vector of SNP identifiers.
}
\author{Robert Scharpf}
\note{
}
\seealso{
  
}

\examples{
data(sample.snpset)
dbSnpId(sample.snpset)[1:5]
}
\keyword{manip}

