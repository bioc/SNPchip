\name{plotCytoband}
\alias{plotCytoband}
\title{Adds cytoband to plot.}

\description{ Adds cytobands to views of chromosome copy number and
  genotype plots.  }

\usage{
plotCytoband(object, ...)
}
\arguments{
  \item{object}{Object of class \code{AnnotatedSnpSet}}
  \item{\dots}{Additional arguments to par such as xlim, cex.axis, xaxs, and main.}
}
\details{
  This function is meant to primarily be called by plotSnp and displayed beneath copy number and genotype estimates.
}

\value{ Only recommended when plotting 5 or fewer chromosomes
(resolution is poor otherwise). It is possible to plot partial views
of a chromosome and its cytoband as well (see examples).  }

\references{}
\author{Robert Scharpf and Jason Ting}
\examples{
data(annSnpset)
data(cytoband)
str(cytoband)
chr1 <- annSnpset[chromosome(annSnpset) == "1", ]

par(mfrow=c(1, 1), mar=c(0, 0, 3, 0), oma=c(5, 0, 3, 0))
plotCytoband(object=chr1, chromosome="1", main="chromosome 1")

plotSnp(chr1, 1, 5, oma=c(6, 4, 0, 3), cex=c(3,3,3),
        cex.axis=1.2, cex.legend=c(1.2, 1.2), addCytoband=TRUE,
        legend.location=c("topleft", "bottomleft"),
        height.cytoband=0.2, width.right=0.2, bty="o", cex.lab=1.5, ncol=1)

##The x-axis begins at the first SNP past the centromere
qarm <- position(chr1) > centromere(chromosomeAnnotation(annSnpset), "1")[[2]] & !is.na(position(chr1))
plotSnp(chr1[qarm, ], 1, 5, oma=c(6, 4, 0, 3), cex=rep(1,3), pch=rep(20,3),
        cex.axis=1.2, cex.legend=c(1.2, 1.2), addCytoband=TRUE,
        legend.location=c("topleft", "bottomleft"),
        height.cytoband=0.2, width.right=0.2, bty="o", cex.lab=1.5, ncol=1, xaxs="r")

##Not equivalent to previous plot because the x-axis begins at the end
##of the centromere (note that there are no SNPs close to the centromere).
qarm <- c(centromere(chromosomeAnnotation(annSnpset), "1")[[2]],
          chromosomeSize(chromosomeAnnotation(annSnpset), "1")[[1]])
plotSnp(chr1, 1, 5, oma=c(6, 4, 0, 3), cex=c(1, 1, 1), pch=rep(20, 3),
        cex.axis=1.2, cex.legend=c(1.2, 1.2), addCytoband=TRUE,
        legend.location=c("topleft", "bottomleft"),
        height.cytoband=0.2, width.right=0.2, bty="o", cex.lab=1.5, ncol=1, adj=0,
        xlim=qarm)
}
\seealso{\code{\link[SNPchip]{plotSnp}}, \code{\link[SNPchip]{AnnotatedSnpSet}}}
\keyword{aplot}

