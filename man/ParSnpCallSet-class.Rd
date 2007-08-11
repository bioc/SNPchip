\name{ParSnpCallSet-class}
\docType{class}
\alias{ParSnpCallSet-class}
\alias{initialize,ParSnpCallSet-method}

\title{Graphical parameters for SnpCopyNumberSet}
\description{contains \code{layout} and \code{par} arguments for \code{SnpCallSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParSnpCallSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}  }
  }
}
\section{Extends}{
  Class \code{\linkS4class{ParESet}}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParSnpCallSet")}: }
  }
}
\references{}
\author{}
\note{}
\seealso{
  \code{\link[oligo]{SnpCopyNumberSet-class}}, \code{\link{par}}, \code{\link[graphics]{layout}}
}
\examples{
showClass("ParSnpCallSet")
}
\keyword{classes}
