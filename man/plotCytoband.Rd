\name{plotCytoband}
\alias{plotCytoband}
\title{Plot the cytobands for a chromosome}
\description{
}
\usage{
plotCytoband(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{AnnotatedSnpSet} object }
  \item{\dots}{Additional arguments}
}
\details{
}
\value{
}
\author{Jason Ting}
\note{
}
\seealso{See Also \code{\link{plotChromosome}}}
\examples{
data(annSnpset)
data(cytoband)
chr1 <- annSnpset[chromosome(annSnpset) == "chr1", 5]
layout(matrix(1:2, nr=2, nc=1), heights=c(1, .05))
par(oma=c(4, 4, 0, 3))
plotChromosome(chr1, mar=c(3, 0, 1, 0), colCentromere="bisque", 
               cexAA=3, cexAB=3, cexNC=3,
               panel.xaxis=TRUE, 
               xlab="Mb", cex.axis=0.7, cex.legend=0.7)
par(las=3)
mtext("copy number", 2, line=2, outer=TRUE, cex=0.7)
par(las=1)
plotCytoband(chr1, cytoband, cex.axis=0.5)
}
\keyword{methods}

