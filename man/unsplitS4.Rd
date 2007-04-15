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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{value}{Object of class \code{AnnotatedSnpSet}}
  \item{featureData}{Object of class \code{AnnotatedDataFrame}}
}
\details{
}
\value{
  \item{\code{AnnotatedSnpSet}}{An object of class \code{AnnoatedSnpSet}}
  ...
}
\references{}
\author{Robert Scharpf}
\note{
}
\seealso{\code{\link[base]{split}}, \code{\link[SNPchip]{AnnotatedSnpSet}}}
\examples{
data(annSnpset)
x <- split(annSnpset, chromosome(annSnpset))
class(x)
length(x)
xx <- unsplitS4(x, featureData(annSnpset))
class(xx)
}
\keyword{manip}

