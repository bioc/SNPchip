\name{AnnotatedSnpSetList-class}
\docType{class}
\alias{AnnotatedSnpSetList-class}
\alias{snpSetList,AnnotatedSnpSetList-method}

\title{Class "AnnotatedSnpSetList"}
\description{Contains a list of AnnotatedSnpSet.  Each element in the list is an AnnotatedSnpSet of one chromosome}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("AnnotatedSnpSetList", snpSetList)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpSetList}}{Object of class \code{list} }
  }
}
\section{Methods}{
  \describe{
    \item{snpSetList}{\code{signature(object = "AnnotatedSnpSetList")}:
      accessor for the class}
  }
}
\references{ }
\author{}
\note{}

\seealso{
}
\examples{
data(hapmap)
df <- data.frame(rep(0, dim(hapmap$calls)[2]), row.names = colnames(hapmap$calls))
colnames(df) <- c("normal")
varMetadata <- data.frame("normal Refset", row.names = "normal")
colnames(varMetadata) <- "labelDescription"
ad <- new("AnnotatedDataFrame", data = df, varMetadata = varMetadata)
snpset <- new("oligoSnpSet", phenoData = ad,
              calls = hapmap$calls,
              callsConfidence = hapmap$callsConfidence,
              cnConfidence = hapmap$callsConfidence, 
              copyNumber = hapmap$copyNumber,
              annotation = "mapping100k")
data(annSnpset)
obj <- as(annSnpset, "AnnotatedSnpSetList")
snpSetList(obj)[1]
snpSetList(obj)[2]
}
\keyword{classes}
