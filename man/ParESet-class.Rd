\name{ParESet-class}
\docType{class}
\alias{ParESet-class}
\alias{$,ParESet-method}
\alias{$<-,ParESet-method}
\alias{allPlots,ParESet-method}
\alias{initialize,ParESet-method}
\alias{snpPar,ParESet-method}
\alias{snpPar<-,ParESet-method}

\title{Class for organizing graphical parameters}
\description{Arguments for plotting high throughput SNP data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParESet", layout, col.axis, cex.main, cex.axis, cex, cex.lab, pch, col, bg, xaxs, xaxt, yaxs, yaxt, lab, adj, bty, ann, useLayout, mar, oma, las, log, ylab, side.ylab, outer.ylab, line.ylab, cex.ylab, xlab, outer.xlab, side.xlab, cex.xlab, line.xlab, outer.axis, line.axis, main, col.centromere, border.centromere, xlim, ylim, one.ylim, add.cytoband, outer.cytoband, outer.cytoband.axis, label.cytoband, use.chromosome.size, label.chromosome, line.label.chromosome, xaxis.side, alternate.xaxis.side, ...)}.
}
\section{Slots}{
  
  \describe{
    \item{\code{snpPar}:}{list containing arguments to \code{par} and \code{layout}}
  }
}
\section{Methods}{
  \describe{
%    \item{$}{\code{signature(x = "ParESet")}  }
%    \item{$<-}{\code{signature(x = "ParESet")}  }
    \item{allPlots}{\code{signature(object = "ParESet")}  }
    \item{initialize}{\code{signature(.Object = "ParESet")}  }
    \item{snpPar}{\code{signature(object = "ParESet")}  }
    \item{snpPar<-}{\code{signature(object = "ParESet")} }
  }
}
\references{}
\author{}
\note{}

\seealso{
  \code{\link{par}}, \code{\link[graphics]{layout}}
}
\examples{
showClass("ParESet")
}
\keyword{classes}
