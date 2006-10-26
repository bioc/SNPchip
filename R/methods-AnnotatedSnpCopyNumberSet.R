setMethod("chromosomeAnnotation", "AnnotatedSnpCopyNumberSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpCopyNumberSet", "data.frame"), function(object, value){
  object@chromosomeAnnotation <- value
  object
})
setMethod("initialize", "AnnotatedSnpCopyNumberSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   copyNumber = new("matrix"),
                   cnConfidence = new("matrix"),
                   chromosomeAnnotation = data.frame())
          {
            .Object@assayData <- assayDataNew(copyNumber = copyNumber, cnConfidence = cnConfidence)
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object
          })
#setMethod("alleleA", "AnnotatedSnpSet", function(object) alleleA(featureData(object)))
#setMethod("alleleB", "AnnotatedSnpSet", function(object) alleleB(featureData(object)))
#setMethod("chromosome", "AnnotatedSnpSet", function(object) chromosome(featureData(object)))
#setMethod("dbSnpId", "AnnotatedSnpSet", function(object) dbSnpId(featureData(object)))
#setMethod("enzyme", "AnnotatedSnpSet", function(object) enzyme(featureData(object)))
#setMethod("position", "AnnotatedSnpSet", function(object) position(featureData(object)))
#setMethod("probeSetId", "AnnotatedSnpSet", function(object) probeSetId(featureData(object)))

