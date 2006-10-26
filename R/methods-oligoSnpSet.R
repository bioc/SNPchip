setAs("oligoSnpSet", "AnnotatedSnpSet", function(from){
  tmp <- new("AnnotatedSnpSet",
             assayData = assayData(from),
             phenoData = phenoData(from),
             featureData = featureData(from),
             experimentData = experimentData(from),
             annotation = annotation(from),
             chromosomeAnnotation = data.frame())
})
setMethod("alleleA", "oligoSnpSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "oligoSnpSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "oligoSnpSet", function(object) chromosome(featureData(object)))
setMethod("dbSnpId", "oligoSnpSet", function(object) dbSnpId(featureData(object)))
setMethod("enzyme", "oligoSnpSet", function(object) enzyme(featureData(object)))
setMethod("position", "oligoSnpSet", function(object) position(featureData(object)))
setMethod("probeSetId", "oligoSnpSet", function(object) probeSetId(featureData(object)))

setMethod("initialize", "oligoSnpSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   featureData = new("AnnotatedDataFrame"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix"),
                   copyNumber = new("matrix"),
                   cnConfidence = new("matrix")){
            .Object@assayData <- assayDataNew(calls = calls,
                                              callsConfidence = callsConfidence,
                                              copyNumber = copyNumber,
                                              cnConfidence = cnConfidence)
            .Object@phenoData <- phenoData
            .Object@experimentData <- experimentData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object
          })
