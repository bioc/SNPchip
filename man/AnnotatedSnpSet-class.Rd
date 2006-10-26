\name{AnnotatedSnpSet-class}
\docType{class}
\alias{AnnotatedSnpSet-class}

\title{Class "AnnotatedSnpSet" }
\description{Container for genotype and copy number data from SNP chips,
as well as chromosome-level features}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AnnotatedSnpSet", phenoData, assayData, experimentData, featureData, annotation, chromosomeAnnotation)}.
}
\section{Slots}{
	 \describe{
    \item{\code{chromosomeAnnotation}:}{Object of class
      \code{data.frame} containing chromosome-level annotation }
    \item{\code{assayData}:}{Object of class \code{AssayData}
      containing copy number and genotype call matrices }
    \item{\code{phenoData}:}{Object of class \code{AnnotatedDataFrame}  }
    \item{\code{featureData}:}{Object of class
      \code{AnnotatedDataFrame}.  Contains annotation for the SNP's,
      such as chromosome, physical position, and alleles. }
    \item{\code{experimentData}:}{Object of class \code{MIAME}  }
    \item{\code{annotation}:}{Object of class \code{character}  }
    \item{\code{.__classVersion__}:}{Object of class \code{Versions}  }
  }
}
\section{Extends}{
Class \code{AnnotatedSnpCallSet}, directly.
Class \code{AnnotatedSnpCopyNumberSet}, directly.
Class \code{SnpCallSet}, by class AnnotatedSnpCallSet, distance 2.
Class \code{eSet}, by class AnnotatedSnpCallSet, distance 3.
Class \code{VersionedBiobase}, by class AnnotatedSnpCallSet, distance 4.
Class \code{Versioned}, by class AnnotatedSnpCallSet, distance 5.
Class \code{SnpCopyNumberSet}, by class AnnotatedSnpCopyNumberSet, distance 2.
Class \code{eSet}, by class AnnotatedSnpCopyNumberSet, distance 3.
Class \code{VersionedBiobase}, by class AnnotatedSnpCopyNumberSet, distance 4.
Class \code{Versioned}, by class AnnotatedSnpCopyNumberSet, distance 5.
}
\section{Methods}{
  \describe{
    \item{alleleB}{\code{signature(object="AnnotatedSnpSet")}: return
      the B allele for one or more SNPs} 
  }
}
\references{ }
\author{Rob Scharpf}
\note{}
\seealso{
}
\examples{
}
\keyword{classes}
