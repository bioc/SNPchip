\name{plotCytoband}
\alias{plotCytoband}
\title{Plots cytoband for one chromosome}

\description{ Adds cytobands to views of chromosome copy number and
  genotype plots.  }

\usage{
plotCytoband(chromosome, cytoband, cytoband.ycoords, xlim, ylim=c(0, 2),
xaxs="r", new=TRUE,
label.cytoband=TRUE, label.y=NULL, srt, cex.axis=1, outer=FALSE, taper=0.15, ...)
}
\arguments{
  \item{chromosome}{character string: which chromosome to draw the
    cytoband}
  \item{cytoband}{data.frame containing cytoband information}
  \item{cytoband.ycoords}{numeric: y coordinates }
  \item{xlim}{x-axis limits}
  \item{ylim}{y-axis limits}
  \item{xaxs}{see par}
  \item{new}{logical: if TRUE, plots the cytoband on a new device}
  \item{label.cytoband}{logical: if TRUE, labels the cytobands}
  \item{label.y}{numeric: height (y-coordinate) for cytoband label}
  \item{srt}{string rotation for cytoband labels.  See \code{par}}
  \item{cex.axis}{size of cytoband labels. See par}
  \item{outer}{logical: whether to draw the labels in the outer
    margins. See par}
  \item{taper}{tapering for the ends of the cytoband}
  \item{\dots}{additional arguments to plot}
}

\author{Robert Scharpf and Jason Ting}
\examples{
plotCytoband("1")
}
\seealso{\code{\link{plotSnp}}, \code{\link{cytoband}}}
\keyword{aplot}

