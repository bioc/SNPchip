\name{xypanel}
\alias{xypanel}
\title{A panel function for plotting copy number versus physical position}
\description{A panel function for \code{xyplot} for plotting copy number versus physical position.}
\usage{
xypanel(x, y, gt, is.snp, range, col.hom = "grey20", fill.hom =
"lightblue", col.het = "grey20", fill.het = "salmon", col.np = "grey20",
fill.np = "grey60", show.state=TRUE, state.cex=1,  col.state="blue", ..., subscripts)
}
\arguments{
  \item{x}{Physical position in megabases.}
  \item{y}{Copy number estimates.}
  \item{gt}{Genotype calls.}
  \item{is.snp}{Logical. Whether the marker is polymorphic.}
  \item{range}{A \code{RangedData} or \code{IRanges} object. Note that we
    expect the units returned by \code{start} and \code{end} to be basepairs.}
  \item{col.hom}{A specification for the color of plotting symbols
    for homozygous genotypes.}
  \item{fill.hom}{A specification for the fill color of plotting symbols
    for homozygous genotypes.
  }
  \item{col.het}{A specification for the color of plotting symbols
    for heterozygous genotypes.}
  \item{fill.het}{A specification for the fill color of plotting symbols
    for heterozygous genotypes.
  }
  \item{col.np}{A specification for the color of plotting symbols
    for nonpolymorphic markers.}
  \item{fill.np}{A specification for the fill color of plotting symbols
    for nonpolymorphic genotypes.
  }
  \item{show.state}{Logical.  Whether to display the predicted state in
    each panel.}
  \item{state.cex}{Numeric.  \code{cex} for state label.  Ignored if
    \code{show.state} is \code{FALSE}.}
  \item{col.state}{Character. color for state label. Ignored if
    \code{show.state} is \code{FALSE}.}
  \item{\dots}{
    Additional arguments passed to lattice functions \code{xyplot},
    \code{lpoints}, and \code{lrect}.
  }
  \item{subscripts}{
    See the panel functions in lattice for more information.
  }
}
\details{

  The order of plotting is (1) nonpolymorphic markers, (2), homozygous
  SNPs, and (3) heterozygous SNPs. Stretches of homozygosity should
  appear as blue using the default color scheme.

}

\note{
  To make the drawing of the \code{range} object border invisible, one can use border="white".
}
\author{
R. Scharpf
}

\seealso{

  \code{\link{xyplot}}

}
\examples{
\dontrun{
	if(require("crlmm") && require("VanillaICE") && require("IRanges")){
    library(lattice)
		library(oligoClasses)
    data(oligoSetExample, package="oligoClasses")
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
		xyplot(cn ~ x | range, data=oligoSet, range=g[2], ##range=fit2[1:10, ],
		       frame=2e6,
		       panel=xypanel, cex=0.3, pch=21, border="blue",
		       scales=list(x="free"),
		       col.hom="lightblue", col.het="salmon", col.np="grey60",
		       fill.np="grey60",
		       xlab="Mb")
		## if xyplot method is masked by lattice, do
		##xyplot <- VanillaICE:::xyplot
	}
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{color}

