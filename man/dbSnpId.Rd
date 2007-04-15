\name{dbSnpId}
\alias{dbSnpId}
\alias{dbSnpId,SnpCallSet-method}
\alias{dbSnpId,SnpCopyNumberSet-method}
\alias{dbSnpId,AnnotatedDataFrame-method}
\alias{dbSnpId,oligoSnpSet-method}
\title{The dbSNP identifier. }
\description{
  Unique identifier for each SNP assigned by NIH
}
\usage{
dbSnpId(object)
}
\arguments{
  \item{object}{Object possibly derived from \code{AnnotatedSnpSet} or
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
\seealso{\code{\link[SNPchip]{AnnotatedSnpSet}},
  \code{\link[Biobase]{featureData}}, \code{\link[Biobase]{AnnotatedDataFrame}}
}
\examples{
data(annSnpset)
dbSnpId(annSnpset)[1:5]
}
\keyword{manip}

