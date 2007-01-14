\name{chromosome-methods}
\docType{methods}
\alias{chromosome-methods}
\alias{chromosome,AnnotatedDataFrame-method}
\alias{chromosome,SnpCallSet-method}
\alias{chromosome,SnpCopyNumberSet-method}
\alias{chromosome,oligoSnpSet-method}
\title{Accessor for SNP annotation}
\description{
accessor for chromosome SNP annotation
}
\section{Methods}{
  \describe{
    
  \item{object = "AnnotatedDataFrame"}{returns the a character string
    for the chromosome number of each SNP}

  \item{object = "SnpCallSet"}{returns the chromosome (character string)
    for each SNP in object of class SnpCallSet }

  \item{object = "SnpCopyNumberSet"}{returns the chromosome (character string)
    for each SNP in object of class SnpCopyNumberSet   }

  \item{object = "oligoSnpSet"}{returns the chromosome (character string)
    for each SNP in object of class oligoSnpSet}

  \item{object = "AnnotatedSnpSet"}{returns the chromosome (character string)
    for each SNP in object of class AnnotatedSnpSet}
}}
\keyword{methods}

