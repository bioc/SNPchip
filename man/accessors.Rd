\name{Accessors}
\alias{chromosome}
\alias{probeSetId}
\alias{dbSnpId}
\title{Accessors for AnnotatedSnpSet objects}
\description{Set of accessors for AnnotatedDataFrame objects.}
\usage{
chromosome(object)
probeSetId(object)
dbSnpId(object)
}
\arguments{
  \item{object}{Object derived from class \code{AnnotatedDataFrame}}
}
\value{
 \code{chromosome} returns the vector of chromosomes for each SNP in the \code{AnnotatedDataFrame} object;
 \code{dbSnpId} returns the dbSNP identifier;
 \code{probeSetId} returns the Affymetrix identifier.
}
\examples{
## Do not run
}
\keyword{manip}% at least one, from doc/KEYWORDS
