\name{getSnpAnnotation}
\alias{getSnpAnnotation}
\alias{getSnpAnnotation,eSet-method}
\title{Retrieves SNP-level annotation}
\description{
  SNP-level annotation currently includes the db SNP identifier,
  chromosome, physical position, strand, and fragment length.
}
\usage{
getSnpAnnotation(object, ...)
}
\arguments{
  \item{object}{An object derived from \code{eSet}.  }
  \item{\dots}{Additional arguments}
}

\details{ Currently, this function supports only Affymetrix 100k and
  500k SNP chips.  The \code{featureNames} in the \code{eSet} object
  must be valid Affymetrix identifiers.  Two options exist for adding
  SNP-level annotation to an object:

  1. The default method is to use the annotation packages maintained by
     Bioconductor.   This is the preferred method as the annotation
     packages are continuously updated.  Using the default method
     requires that the annotation slot in the \code{eSet} object be
     exactly identical to 1 of the 4 annotation packages currently
     provided:

     a. pd.mapping50k.hind240
     b. pd.mapping50k.xba240
     c. pd.mapping250k.nsp
     d. pd.mapping250k.sty

     Note that if the R package oligo was used to summarize genotype
     calls in the \code{eSet} object, the annotation slot is already
     provided.

  2. Manually download the annotation provided by NetAffx and convert
     the file to an object of class \code{AnnotatedDataFrame}.  The SNPs
     in the \code{AnnotatedDataFrame} must be in the same order as the
     SNPs in the object derivied from \code{eSet}.
}
  
\value{ An object of class \code{AnnotatedDataFrame} } \references{}
  \author{Robert Scharpf} \seealso{\code{\link[oligo]{oligoSnpSet}},
  \code{\link[Biobase]{AnnotatedDataFrame}},
  \code{\link[Biobase]{annotation}}}

\examples{

##Default method for adding SNP-level annotation

  data(annSnpset)
  obj <- new("AnnotatedSnpCopyNumberSet",
             copyNumber=copyNumber(annSnpset)[1:50,],
             cnConfidence=cnConfidence(annSnpset)[1:50,],
             phenoData=phenoData(annSnpset), annotation=annotation(annSnpset))
  obj
\dontrun{
##Requires pd.mapping50k.xba240 annotation package
  snp.annotation <- getSnpAnnotation(obj)
  snp.annotation
  class(snp.annotation)
  identical(featureNames(obj), featureNames(snp.annotation))
  featureData(obj) <- snp.annotation
  validObject(obj)
  varLabels(featureData(obj))
}

##Using NetAffx annotation (already saved as a data.frame)
\dontrun{
  try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kXba240.rda")))
  reorderRows <- match(featureNames(obj), as.character(mapping50kXba240$Probe.Set.ID))
  netaffx.annotation <- mapping50kXba240[reorderRows, ]
  fD <- as(netaffx.annotation, "AnnotatedDataFrame")
  featureData(obj) <- fD
  varLabels(featureData(obj))
  validObject(obj)
}
}

\keyword{methods}
\keyword{manip}
