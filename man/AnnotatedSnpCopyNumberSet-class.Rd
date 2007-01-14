\name{AnnotatedSnpCopyNumberSet-class}
\docType{class}
\alias{AnnotatedSnpCopyNumberSet-class}
\alias{plotCytoband,AnnotatedSnpCopyNumberSet-method}
\alias{plotChromosome,AnnotatedSnpCopyNumberSet-method}
\alias{plotSnp,AnnotatedSnpCopyNumberSet-method}
\alias{snpSetList,AnnotatedSnpCopyNumberSet-class}

\title{Class "AnnotatedSnpCopyNumberSet"}
\description{Container for copy number data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AnnotatedSnpCopyNumberSet", phenoData, experimentData, annotation, copyNumber, cnConfidence, snpLocation, chromosomeAnnotation, chromosome)}.
}
\section{Slots}{
  \describe{
    \item{\code{chromosomeAnnotation}:}{Object of class \code{data.frame}  }
    \item{\code{assayData}:}{Object of class \code{AssayData}  }
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame} }
    \item{\code{experimentData}:}{Object of class \code{MIAME}  }
    \item{\code{featureData}:}{Object of class  \code{AnnotatedDataFrame}.  Contains annotation for the SNP's,
      such as chromosome, physical position, and alleles. }    
    \item{\code{annotation}:}{Object of class \code{character}  }
    \item{\code{.__classVersion__}:}{Object of class \code{Versions} }
  }
}
\section{Extends}{
Class \code{SnpCopyNumberSet}, directly.
Class \code{eSet}, by class \code{SnpCopyNumberSet}.
Class \code{VersionedBiobase}, by class \code{SnpCopyNumberSet}.
Class \code{Versioned}, by class \code{SnpCopyNumberSet}.
}
\section{Methods}{
  \describe{
    Methods for AnnotatedSnpCopyNumberSet
  }
}
\references{}
\author{}
\note{}
\seealso{
}
\examples{
}
\keyword{classes}
