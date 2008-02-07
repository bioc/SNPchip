\name{chromosome2numeric}
\alias{chromosome2numeric}
\title{Converts chromosome to numeric}
\description{
  Provides a standard way of converting chromosome annotation stored in
  \code{SnpLevelSet} objects to the numeric class.
}
\usage{
chromosome2numeric(chromosome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chromosome}{ an object inherited from \code{SnpLevelSet}}
}
\details{
  This is useful when sorting SNPs in an object by chromosome and
  physical position -- ensures that the sorting is done in the same way
  for different objects.

  The function currently provides the following conversions:

  "X" -> 23
  "XY" -> 24
  "Y" -> 25
  "M" -> 26
  }

\value{
  \code{numeric}
}
\author{R. Scharpf}
\keyword{manip}

