\name{AnnotatedSnpSet-class}
\docType{class}
\alias{AnnotatedSnpSet-class}
\alias{alleleA,AnnotatedSnpSet-method}
\alias{alleleB,AnnotatedSnpSet-method}
\alias{chromosome,AnnotatedSnpSet-method}
\alias{enzyme,AnnotatedSnpSet-method}
\alias{plotCytoband,AnnotatedSnpSet-method}
\alias{plotSnp,AnnotatedSnpSet-method}
\alias{position,AnnotatedSnpSet-method}
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
    \item{alleleA}{\code{signature(object="AnnotatedSnpSet")}: return
      the A allele for the specified SNPs in object}
    
    \item{alleleB}{\code{signature(object="AnnotatedSnpSet")}: return
      the B allele for the specified SNPs in object}

    \item{chromosome}{\code{signature(object="AnnotatedSnpSet")}:
      accessor for the chromosome}

    \item{plotCytoband}{\code{signature(object="AnnotatedSnpSet")}:
    Plots the cytoband for a single chromosome.}
  
    \item{plotSnp}{\code{signature(object="AnnotatedSnpSet")}: Makes a
      scatterplot of copy number (vertical axis) versus physical
      position (horizontal axis) for each chromosome and sample in the
      AnnotatedSnpSet object.  Points are color-coded by the genotype
      call.}

    \item{position}{\code{signature(object="AnnotatedSnpSet")}: Returns
      physical position of SNP along the chromosome in numbers of
      basepairs from the 5' end (an integer)}
  }
}
\references{ }
\author{Rob Scharpf}
\note{}
\seealso{
}
\examples{
data(annSnpset)
annSnpset
##for formatting the plot layout, see the vignette
plotSnp(annSnpset, 1:10, 1:2)
}
\keyword{classes}
