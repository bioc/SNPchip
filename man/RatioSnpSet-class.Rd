\name{RatioSnpSet-class}
\docType{class}
\alias{RatioSnpSet-class}
\title{Class for high throughput SNP-level data.  The assay data
  contains ratios of copy number.}
\description{Container for high throughput SNP data for platforms that
  provide ratios as SNP-level summaries of copy number}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RatioSnpSet", assayData, phenoData, featureData, experimentData, annotation, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"} }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpLevelSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "SnpLevelSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpLevelSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "SnpLevelSet", distance 4.
}
\section{Methods}{
No methods defined with class "RatioSnpSet" in the signature.
}
\author{RS}
\note{Might be useful for bead-studio processed data}
\examples{
showClass("RatioSnpSet")
}
\keyword{classes}
