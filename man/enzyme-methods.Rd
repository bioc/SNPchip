\name{enzyme-methods}
\docType{methods}
\alias{enzyme}
\alias{enzyme-methods}
\alias{enzyme,AnnotatedDataFrame-method}
\alias{enzyme,SnpCopyNumberSet-method}
\alias{enzyme,SnpCallSet-method}
\alias{enzyme,oligoSnpSet-method}
\alias{enzyme,AnnotatedSnpSet-method}
\title{Enzyme for DNA digestion}
\description{
Accessor to retrieve which enzyme was used to produce the probe
fragments that comprise the probe set of a SNP
}
\section{Methods}{
  \describe{
    
    \item{object = "AnnotatedDataFrame"}{enzyme used to digest the DNA}
  
    \item{object = "SnpCopyNumberSet"}{ enzyme used to digest the DNA }
    
    \item{object = "SnpCallSet"}{ enzyme used to digest the DNA }

    \item{object = "oligoSnpSet"}{ enzyme used to digest the DNA }

    \item{object = "AnnotatedSnpSet"}{ enzyme used to digest the DNA }
}}
\keyword{methods}

