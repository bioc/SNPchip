\name{xyplot}
%\alias{xyplot}
\alias{xyplot2}
%\alias{xyplot,formula,gSet-method}
%\alias{xyplot,formula,BeadStudioSet-method}
\alias{xyplot2,formula,gSet-method}
\alias{xyplot2,formula,SnpSet-method}
\title{Plot copy number and physical position for a set of genomic intervals.}
\description{Plot copy number and physical position given by a
  \code{CNSet} object for a set of genomic intervals stored in a
  \code{RangedDataCVN} object.
}
\usage{
%xyplot(x, data, ...)
xyplot2(x, data, range, frame=50e3L, ...)
}
\arguments{
  \item{x}{A \code{formula}.  Currently, the formula must be one of \code{cn~x},
    \code{cn ~ x | id} or \code{cn ~ x | range} when \code{data} is a
    \code{CNSet}. If \code{data} is a \code{BeadStudioSet}, the formula
    has the form \code{lrr ~ x| range} or \code{baf ~ x | range}.}
  \item{data}{A \code{CNSet}, \code{BeadStudioSet}, or \code{SnpSet} object.}
  \item{\dots}{A \code{RangedDataCNV} object must be passed by the name
    'range'. Arguments for \code{xyplot} are passed to
    \code{xyplot2}. Additional arguments are passed to \code{xypanel} and \code{panel.xyplot}.}
  \item{range}{A \code{RangedDataCNV} object.}
  \item{frame}{The genomic distance (basepairs) to the left and right of the start and stop
    coordinates in the \code{range} object.}
}
\details{


  These functions plot copy number estimates versus physical
  position. The function is particularly useful for multi-panel displays
  in which the copy number estimates for a single range of a
  \code{GRanges} object appears in one panel. The size of the
  multi-panel display depends on the number of ranges (rows) in the
  \code{GRanges} object.

}

\value{

  An object of class \code{trellis}.

}

\author{
R. Scharpf
}

\seealso{

	\code{\link[lattice]{xyplot}}, \code{\link{xypanel}}

  To modify the plot appearance from the default, additional arguments
  can be passed to \code{\link{panel.xyplot}}, \code{\link{lpoints}}, and
  \code{\link{lrect}}.

}
\examples{
  ## simulated data
  library(oligoClasses)
  library(IRanges)
  library(VanillaICE)
  data(oligoSetExample, package="oligoClasses")
  ## The oligoSnpSet class will likely be deprecated and made defunct
  ## in a future release.  Instead, we favor
  ## RangedSummarizedExperiment-derived classes defined in VanillaICE
  oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
  cn <- copyNumber(oligoSet)/100
  cn <- log2((2^cn)/2)
  gt <- calls(oligoSet)[,]
  ## simulate BAFs
  bf <- rep(NA, length(gt))
  u <- runif(length(gt))
  bf[gt==1 & u > 0.5] <- runif(sum(gt==1 & u > 0.5), 0, 0.05)
  bf[gt==1 & u <= 0.5] <- runif(sum(gt==1 & u <= 0.5), 0.95, 1)
  bf[gt==2] <- runif(sum(gt==2), 0.45, 0.55)
  bf[900:1200] <- runif(length(900:1200), 0, 0.03)
  gr <- GRanges(paste0("chr", chromosome(featureData(oligoSet))),
                IRanges(position(oligoSet), width=1))
  cn <- as.matrix(cn)
  bf <- as.matrix(bf)
  dimnames(cn) <- dimnames(bf) <- list(featureNames(oligoSet), sampleNames(oligoSet))
  se <- SnpArrayExperiment(cn=cn,
                           baf=bf,
                           rowRanges=gr,
                           isSnp=rep(TRUE, length(gr)))
  fit <- hmm2(se)
  g <- as(segs(fit), "GRanges")
  ## To visualize each range in it's own panel surrounded by a
  ## frame of 2,000,000 bases:
  ## (here the frames are overlapping, but the method could be
  ## applied more generally to a collection of ranges from
  ## different chromsomes and samples)
  xyplot2(cn~x | range, data=oligoSet,
	       range=g,
	       frame=2e6, panel=xypanel,
	       cex=2,
	       pch=".",
	       col.het="salmon",
	       fill.het="salmon",
	       col.hom="royalblue",
	       fill.hom="royalblue",
	       state.cex=0.5,
	       border="orange", scales=list(x="free"),
	       par.strip.text=list(cex=0.5),
	       xlab="Mb", ylab=expression(log[2]("copy number")))
}
\keyword{dplot}
\keyword{methods}
