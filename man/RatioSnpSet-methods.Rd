\name{RatioSnpSet-methods}
\alias{initialize,RatioSnpSet-method}
\alias{ratio}
\alias{ratio<-}
\alias{ratio,RatioSnpSet-method}
\alias{ratio,SnpLevelSet-method}
\alias{ratio<-,SnpLevelSet,matrix-method}
\alias{copyNumber}
\alias{copyNumber,RatioSnpSet-method}
\title{Accessors for RatioSnpSet objects}
\description{
  Accessors and replacement methods for RatioSnpSet objects.
}
\usage{
%copyNumber(object)
ratio(object)
ratio(object)<-value
}
\arguments{
  \item{object}{\code{RatioSnpSet} object}
  \item{value}{matrix}
}
\value{
	
  Matrix with R rows (SNPs) and C columns (Samples) containing ratios
  of allele intensities

}

\keyword{manip}
