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
  Class \code{\link[oligoClasses:SnpLevelSet-class]{SnpLevelSet}}, directly.
  Class \code{NULLorHmmPredict}, directly.
  Class \code{\link[Biobase:class.eSet]{eSet}},  by class "SnpLevelSet", distance 2.
  Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "SnpLevelSet", distance 3.
  Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "SnpLevelSet", distance 4.
}
\section{Methods}{
  No methods defined with class "HmmPredict" in the signature.
}
\author{RS}
\examples{
showClass("HmmPredict")
}
\keyword{classes}
