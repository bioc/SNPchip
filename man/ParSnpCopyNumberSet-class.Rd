\name{ParSnpCopyNumberSet-class}
\docType{class}
\alias{ParSnpCopyNumberSet-class}
\alias{initialize,ParSnpCopyNumberSet-method}
\title{Graphical parameters for SnpCopyNumberSet}
\description{contains \code{layout} and \code{par} arguments for \code{SnpCopyNumberSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParSnpCopyNumberSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}.  Contains
      arguments for \code{par} and \code{layout} }
  }
}
\section{Extends}{
Class \code{\linkS4class{ParESet}}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParSnpCopyNumberSet")}: ... }
  }
}
\references{}
\author{}
\note{}

\seealso{
  \code{\link[oligo]{SnpCopyNumberSet-class}}, \code{\link{par}}, \code{\link[graphics]{layout}}
}
\examples{
showClass("ParSnpCopyNumberSet")
}
\keyword{classes}
