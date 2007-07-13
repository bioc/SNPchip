\name{AnnotatedSnpCallSet-class}
\docType{class}
\alias{AnnotatedSnpCallSet-class}
\alias{alleleA,SnpCallSet-class}
\alias{alleleB,SnpCallSet-class}
\alias{chromosomeAnnotation,AnnotatedSnpCallSet-method}
\alias{chromosomeAnnotation<-,AnnotatedSnpCallSet,data.frame-method}


\title{Class "AnnotatedSnpCallSet"}
\description{Container for genotype call data accompanied by SNP annotation}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AnnotatedSnpCallSet", phenoData, experimentData, annotation, calls, callsConfidence, snpLocation, chromosomeAnnotation, chromosome)}.
}
\section{Slots}{
  \describe{
    \item{\code{chromosomeAnnotation}:}{Object of class \code{data.frame}}
    \item{\code{assayData}:}{Object of class \code{AssayData}}
    \item{\code{phenoData}:}{Object of class  \code{AnnotatedDataFrame}}
    \item{\code{featureData}:}{Object of class
      \code{AnnotatedDataFrame} contains annotation for the SNP's, such as chromosome, physical position, and alleles. }    
    \item{\code{experimentData}:}{Object of class \code{MIAME}}
    \item{\code{annotation}:}{For 100k chips, "mapping50kHind240",
      "mapping50kXba240" , or "mapping100k".  For 500k chips, mmust be
      "mapping250kNsp" or "mapping250kSty"}
    \item{\code{.__classVersion__}:}{Object of class \code{Versions}}
  }
}
\section{Extends}{
Class \code{SnpCallSet}, directly.
Class \code{eSet}, by class \code{SnpCallSet}.
Class \code{VersionedBiobase}, by class \code{SnpCallSet}.
Class \code{Versioned}, by class \code{SnpCallSet}.
}
\section{Methods}{
  \describe{
  }
}
\author{Robert Scharpf}
\note{}

\seealso{\code{\link{AnnotatedSnpCopyNumberSet-class}},
  \code{\link{AnnotatedSnpSet-class}}
}
\examples{
}
\keyword{classes}
