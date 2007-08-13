\name{unsplitSnpSet}
\alias{unsplitSnpSet}
\alias{unsplitSnpSet-methods}
\alias{unsplitSnpSet,list,AnnotatedDataFrame-method}
\title{Combines a list with elements that are AnnotatedSnpSets into one AnnotatedSnpSet.}
\description{
  After splitting a \code{SnpSet} object by chromosome, for instance,
  this function will reproduce the original object from the list.
}
\usage{
unsplitSnpSet(from, annotatedDataFrame, ...)
}
\arguments{
  \item{from}{A \code{list}. Element should be one of the SNP classes
    extending \code{eSet}}
  \item{annotatedDataFrame}{Object of class \code{AnnotatedDataFrame} that
    contains the \code{featureData} prior to splitting.  See the
    example.}
  \item{\dots}{arguments to initialize an instance of the class}
}
\details{
}
\value{
  an object of class \code{class(from[[1]])}

}
\references{}
\author{Robert Scharpf}
\note{
}
\seealso{
  \code{\link[base]{split}}
  }
\examples{
data(sample.snpset)
x <- split(sample.snpset, chromosome(sample.snpset))
class(x)
length(x)
xx <- unsplitSnpSet(x, featureData(sample.snpset),
                    copyNumber=do.call("rbind", lapply(x, copyNumber)),
                    calls=do.call("rbind", lapply(x, calls)),
                    phenoData=phenoData(x[[1]]),
                    annotation=annotation(x[[1]]),
                    experimentData=experimentData(x[[1]]))
class(xx)
##is FALSE
identical(featureNames(sample.snpset), featureNames(xx))

##but
sample.snpset <- sample.snpset[order(position(sample.snpset)), ]
xx <- xx[order(position(xx)), ]
identical(featureNames(sample.snpset), featureNames(xx))
}
\keyword{manip}

