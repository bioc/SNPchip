\name{fData}
\alias{fData}
\alias{fData,eSet-method}
\alias{fData<-}
\alias{fData<-,eSet,data.frame-method}
\title{Accessor for the data.frame in featureData}

\description{ This generic function accesses the \code{data.frame} for
the feature-level data, such as annotation of SNPs or genes.  fData is
the analog of \code{pData} that accesses the \code{data.frame} for the
\code{phenoData} \code{\link[Biobase]{AnnotatedDataFrame-class}}.
}

\usage{
fData(object)
}
\arguments{
  \item{object}{Object, possibly derived from class \code{eSet}}
}
\details{
}

\value{
  \code{fData} and \code{pData} both return data.frames.  \code{fData}
  will generally return a very large data.frame as the number of
  features is often very large.
}

\references{}
\author{Robert Scharpf}
\note{}
\seealso{
  \code{\link[Biobase]{AnnotatedDataFrame-class}}
}
\examples{
data(annSnpset)
fData(annSnpset)[1:5, ]
}
\keyword{manip}

