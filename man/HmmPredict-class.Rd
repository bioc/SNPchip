\name{HmmPredict-class}
\docType{class}
\alias{HmmPredict-class}
\title{Class for containing output from a hidden Markov fit to SNP-level
  data of genotype calls and copy number estimates}
\description{Container for HMM predictions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HmmPredict", assayData, phenoData, featureData, experimentData, annotation, ...)}.

}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"character"}  }
    \item{\code{breakpoints}:}{Object of class \code{"data.frame"} }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpLevelSet}"}, directly.
Class \code{NULLorHmmPredict}, directly.
Class \code{"\linkS4class{eSet}"}, by class "SnpLevelSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpLevelSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "SnpLevelSet", distance 4.
}
\section{Methods}{
No methods defined with class "HmmPredict" in the signature.
}
\author{RS}
\examples{
showClass("HmmPredict")
}
\keyword{classes}
