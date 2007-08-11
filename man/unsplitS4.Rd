\name{unsplitS4}
\alias{unsplitS4}
\title{Combines a list with elements that are AnnotatedSnpSets into one AnnotatedSnpSet.}
\description{
  After splitting an \code{AnnotatedSnpSet} by chromosome, for instance,
  this function will reproduce the original object from the list.
}
\usage{
unsplitS4(value, featureData)
}
\arguments{
  \item{value}{A \code{list}. Each element in the must be of class \code{oligoSnpSet}}
  \item{featureData}{Object of class \code{AnnotatedDataFrame} that
    contains the \code{featureData} prior to splitting.  See the example.}
}
\details{
}
\value{
  \code{oligoSnpSet}
  ...
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
xx <- unsplitS4(x, featureData(sample.snpset))
class(xx)
}
\keyword{manip}

