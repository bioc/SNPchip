\name{plotCytoband}
\alias{plotCytoband}
\title{Plots cytoband for one chromosome}

\description{ Adds cytobands to views of chromosome copy number and
  genotype plots.  }

\usage{
plotCytoband(chromosome, cytoband, xlim, xaxs="r", new=TRUE,
label.cytoband=TRUE, cex.axis=1, outer=FALSE, ...)
}
\arguments{
  \item{chromosome}{character string: which chromosome to draw the
    cytoband}
  \item{cytoband}{data.frame containing cytoband information}
  \item{xlim}{axis limits}
  \item{xaxs}{see par}
  \item{new}{logical: if TRUE, plots the cytoband on a new device}
  \item{label.cytoband}{logical: if TRUE, labels the cytobands}
  \item{cex.axis}{size of cytoband labels. See par}
  \item{outer}{logical: whether to draw the labels in the outer
    margins. See par}
  \item{\dots}{additional arguments to plot}
}

\details{

  This function is meant to primarily be called by plotSnp and displayed
  beneath copy number and genotype estimates.
  
}

\value{

}

\references{}
\author{Robert Scharpf and Jason Ting}
\examples{
source("~/projects/software/SNPchip/R/functions.R")
plotCytoband("1")
}
\seealso{\code{\link{plotSnp}}, \code{\link{cytoband}},
  \code{\link{par}}, \code{\link{plot}}}
\keyword{aplot}

