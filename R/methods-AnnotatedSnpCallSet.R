setMethod("chromosomeAnnotation", "AnnotatedSnpCallSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpCallSet", "data.frame"), function(object, value){
  object@chromosomeAnnotation <- value
  object
})
setMethod("initialize", "AnnotatedSnpCallSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix"),
                   chromosomeAnnotation = data.frame()){
            .Object@assayData <- assayDataNew(calls = calls, callsConfidence = callsConfidence)
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object
          })
