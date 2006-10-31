\name{plotSnp-methods}
\docType{methods}
\alias{plotSnp-methods}
\alias{plotSnp,AnnotatedSnpSet-method}
\alias{plotSnp}
\title{Methods for function plotSnp}
\description{
Methods for function \code{plotSnp} in Package SNPchip. Produces a
scatterplot of copy number (vertical axis) versus physical position
(horizontal axis) for each chromosome and sample in the AnnotatedSnpSet
object.  Points are color-coded by the genotype call.
}
\section{Methods}{
\describe{
\item{object = "AnnotatedSnpSet"}{an object of class AnnotatedSnpSet}
}}
\examples{
data(annSnpset)
##for formatting the plot layout, see the vignette
plotSnp(chromosomes = 1:10, object = annSnpset, samples = 1:2, height.bottom =
.5)
}
\keyword{methods}


