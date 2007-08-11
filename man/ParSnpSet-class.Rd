\name{ParSnpSet-class}
\docType{class}
\alias{ParSnpSet-class}
\alias{initialize,ParSnpSet-method}

\title{Graphical parameters for oligoSnpSet}
\description{}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParSnpSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}  }
  }
}
\section{Extends}{
Class \code{\linkS4class{ParSnpCopyNumberSet}}, directly.
Class \code{\linkS4class{ParESet}}, by class "ParSnpCopyNumberSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParSnpSet")}: ... }
  }
}
\references{}
\author{}
\note{}

\seealso{
  \code{\link[oligo]{oligoSnpSet-class}}, \code{\link{par}}, \code{\link[graphics]{layout}}
}
\examples{
showClass("ParSnpSet")
}
\keyword{classes}
