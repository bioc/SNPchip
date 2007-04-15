\name{plotSnp}
\alias{plotSnp}
\title{Plots copy number and genotype calls.}
\description{Plots copy number and genotype call versus physical position for
  an arbitrary number of samples and chromosomes.
}
\usage{
plotSnp(object, chromosomes, samples, snpId = c(NA, 1e+06), oma = c(5, 4, 4, 0.5), mar = c(0.5, 0, 0.5, 0.2), width.right = NULL, showLayout = FALSE, plot = TRUE, col = c("royalblue", "red", "royalblue", "green3"), bg = rep("white", 4), cex = c(2, 3, 2, 2), pch = rep(".", 4), col.centromere = c("bisque", NA), bw = FALSE, lwd = 2, col.axis = "brown", cex.axis = 1, cex.main = 1, cex.lab = 1, xlim = NULL, ylim = NULL, log = "", xaxis.side = rep(c(1, 3), length.out = length(chromosomes)), xaxs = "i", xaxt = "n", yaxs = "i", yaxt = "n", xlab = "", ylab = "", lab = c(2, 5, 7), adj = 0, main = "", legend = c(TRUE, TRUE), legend.panel = c(TRUE, FALSE), legend.location = c("topright", "bottomright"), legend.bty = c("n", "n"), legend.col = c("white", "white"), ncol = 2, cex.legend = c(1, 1), digits = 3, pt.cex = cex.legend * 1.5, bty = "n", addCytoband = FALSE, height.cytoband = 0.5, ...)
}
\arguments{
  \item{object}{Object of class \code{AnnotatedSnpCopyNumberSet} or \code{AnnotatedSnpSet}}
  \item{chromosomes}{Vector of character strings for chromosome number.
  Chromosomes will be plotted in the order given.}
  \item{samples}{Vector of integers indicating sample number.}
  \item{snpId}{Not yet implemented.}
  \item{oma}{See \code{par}}
  \item{mar}{See \code{par}}
  \item{width.right}{Width of panel for the legend relative to the width
    of the smallest chromosome.}
  \item{showLayout}{Logical for whether to display the output of
    \code{layout} for plotting}
  \item{plot}{Logical.  If FALSE, only the output from \code{layout} is displayed.}
  \item{col}{A character or numerical vector indicating the color to plot genotype
    call:  AA color, AB color,  BB color, no call color.  For
    \code{AnnotatedSnpCopyNumberSet}, only the first color is used.}
  \item{bg}{A character or numerical vector indicating the background color for the
    plotting symbols.  Possible to specify a different background for
    each genotype call in the same order as \code{col}.}
  \item{cex}{A numerical vector of plotting symbol size in the same
    order as \code{col}.}
  \item{pch}{A character or numerical vector of plotting symbol size in the same
    order as \code{col}.}
  \item{col.centromere}{Vector of length two: color of centromere and
    color for the centromere border. }
  \item{bw}{Converts all colors to grayscale.  }
  \item{lwd}{See \code{par}}
  \item{col.axis}{See \code{par}}
  \item{cex.axis}{See \code{par}}
  \item{cex.main}{See \code{par}}
  \item{cex.lab}{See \code{par}}
  \item{xlim}{See \code{par}}
  \item{ylim}{See \code{par}}
  \item{log}{See \code{par}}
  \item{xaxis.side}{Vector of integers indicating whether to plot the
    x-axis on the top (3) or bottom (1) for each sample.  The default is
  to alternate top/bottom.}
  \item{xaxs}{See \code{par}}
  \item{xaxt}{See \code{par}}
  \item{yaxs}{See \code{par}}
  \item{yaxt}{See \code{par}}
  \item{xlab}{See \code{par}}
  \item{ylab}{See \code{par}}
  \item{lab}{See \code{par}}
  \item{adj}{See \code{par}}
  \item{main}{See \code{par}}
  \item{legend}{Logical vector of length 2.  If both are TRUE, a legend
    for descriptive statistics (legend 1) and plotting symbols (legend
    2) are generated.}
  \item{legend.panel}{Logical vector of length 2.  Whether to draw the
    legends on a separate panel.}
  \item{legend.location}{Character vector of length 2.  Where to draw
    the legends. See \code{legend}.}
  \item{legend.bty}{Character vector of length 2.  See the argument
    \code{bty} in function \code{legend}.}
  \item{legend.col}{Character vector of length 2. See the argument
    \code{bty} in function \code{legend}. }
  \item{ncol}{Number of columns in legend for plotting symbols.  See \code{legend}.}
  \item{cex.legend}{Same as argument \code{cex} in function \code{legend}.}
  \item{digits}{Rounds values in legend to the number of digits specified.}
  \item{pt.cex}{See \code{legend}}
  \item{bty}{See \code{par}}
  \item{addCytoband}{Logical indicating whether to add a cytoband.}
  \item{height.cytoband}{Height of cytoband relative to height of samples.}
  \item{\dots}{ Additional arguments}
}
\details{
}
\value{
}
\references{}
\author{Robert Scharpf}
\note{
}
\seealso{\code{\link{par}}, \code{\link{plot}}, \code{\link{legend}}}
}
\examples{
data(annSnpset)
plotSnp(annSnpset, chromosomes=c(1:22, "X"), 
        samples=1:4, width.right=10, cex.axis=1, lab=c(3, 3, 5),
        cex.lab=1.2, cex.legend=c(1, 1.2))
}
\keyword{hplot}

