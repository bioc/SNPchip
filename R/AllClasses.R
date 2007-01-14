setClass("AnnotatedSnpCallSet", representation(chromosomeAnnotation="data.frame"), contains="SnpCallSet")
setClass("AnnotatedSnpCopyNumberSet", representation(chromosomeAnnotation="data.frame"), contains="SnpCopyNumberSet")
setClass("AnnotatedSnpSet",  contains = c("AnnotatedSnpCallSet", "AnnotatedSnpCopyNumberSet"))
setClass("AnnotatedSnpSetList", representation(snpSetList="list"))
setClass("AnnotatedSnpCopyNumberSetList", representation(snpSetList="list"))
setClass("AnnotatedSnpCallSetList", representation(snpSetList="list"))

setValidity("AnnotatedSnpSet", function(object){
  assayDataValidMembers(assayData(object), c("calls", "copyNumber", "callsConfidence", "cnConfidence"))
})
setValidity("AnnotatedSnpCallSet", function(object){
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence"))
})
setValidity("AnnotatedSnpCopyNumberSet", function(object){
  assayDataValidMembers(assayData(object), c("copyNumber", "cnConfidence"))
})




