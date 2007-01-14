\name{AnnotatedSnpCopyNumberSetList-class}
\docType{class}
\alias{AnnotatedSnpCopyNumberSetList-class}
\alias{coerce,AnnotatedSnpCopyNumberSet,AnnotatedSnpCopyNumberSetList-method}
\alias{coerce,AnnotatedSnpCopyNumberSetList,AnnotatedSnpCopyNumberSet-method}
\alias{initialize,AnnotatedSnpCopyNumberSetList-method}
\alias{snpSetList,AnnotatedSnpCopyNumberSetList-method}

\title{A list with elements of class AnnotatedSnpCopyNumberSet}
\description{Each element of the list is typically an
  AnnotatedSnpCopyNumberSet of one chromosome}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{as(object,
    "AnnotatedSnpCopyNumberSetList")}, where object is of class \code{AnnotatedSnpCopyNumberSet}.
}
\section{Slots}{
	 \describe{
    \item{\code{snpSetList}:}{Object of class \code{"list"}}
  }
}
\section{Methods}{
  \describe{
    
    \item{snpSetList}{\code{signature(object="AnnotatedSnpCopyNumberSetList")}:
      extracts list of elements \code{AnnotatedSnpCopyNumberSet}}

  }
    
}
\author{Robert Scharpf}
\note{}
\seealso{
	\code{\linkS4class{AnnotatedSnpCopyNumberSet}} 
}
\examples{
}
\keyword{classes}
