\name{RatioSnpSet-methods}
\alias{RatioSnpSet-methods}
\alias{initialize,RatioSnpSet-method}
%\alias{cnConfidence}
%\alias{cnConfidence<-}
%\alias{copyNumber}
%\alias{copyNumber<-}
\alias{ratio}
\alias{ratio<-}
\alias{ratio,RatioSnpSet-method}
\alias{ratio,SnpLevelSet-method}
\alias{ratio<-,SnpLevelSet,matrix-method}
\alias{cnConfidence,RatioSnpSet-method}
\alias{cnConfidence<-,RatioSnpSet,matrix-method}
\alias{copyNumber,RatioSnpSet-method}
\alias{copyNumber<-,RatioSnpSet,matrix-method}
\alias{ratioConfidence}
\alias{ratioConfidence<-}
\alias{ratioConfidence,RatioSnpSet-method}
\alias{ratioConfidence<-,RatioSnpSet,matrix-method}
\title{Accessors for RatioSnpSet objects}
\description{
  Accessors and replacement methods for RatioSnpSet objects.
}
\usage{
	ratio(object)
	ratio(object)<-value
	ratioConfidence(object)
	ratioConfidence(object) <- value
%	cnConfidence(object)
%	cnConfidence(object) <- value
%	copyNumber(object)
%	copyNumber(object) <- value
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
