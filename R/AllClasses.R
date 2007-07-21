setClass("AnnotatedSnpCallSet",
         representation(chromosomeAnnotation="data.frame"),
         prototype=list(
           chromosomeAnnotation=data.frame()),
         contains="SnpCallSet")

setClass("AnnotatedSnpCopyNumberSet",
         representation(chromosomeAnnotation="data.frame"),
         prototype=list(
           chromosomeAnnotation=data.frame()),         
         contains="SnpCopyNumberSet")

setClass("AnnotatedSnpSet",
         contains = c("AnnotatedSnpCallSet", "AnnotatedSnpCopyNumberSet"))

setValidity("AnnotatedSnpSet", function(object){
  assayDataValidMembers(assayData(object), c("calls", "copyNumber", "callsConfidence", "cnConfidence"))
})

setValidity("AnnotatedSnpCallSet", function(object){
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence"))
})

setValidity("AnnotatedSnpCopyNumberSet", function(object){
  assayDataValidMembers(assayData(object), c("copyNumber", "cnConfidence"))
})

setAs("AnnotatedSnpSet", "AnnotatedSnpCopyNumberSet",
      function(from){
        object <- new("AnnotatedSnpCopyNumberSet",
                      copyNumber=copyNumber(from),
                      cnConfidence=cnConfidence(from),
                      phenoData=phenoData(from),
                      featureData=featureData(from),
                      annotation=annotation(from),
                      chromosomeAnnotation=chromosomeAnnotation(from))
        object
      })




