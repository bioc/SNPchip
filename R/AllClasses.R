setClass("ParESet", representation(snpPar="list"))
setClass("ParSnpCopyNumberSet", contains="ParESet")
setClass("ParSnpCallSet", contains="ParESet")
setClass("ParSnpSet", contains="ParSnpCopyNumberSet")

###########################################################################
##DEPRECATED CLASSES
###########################################################################
setClass("AnnotatedSnpSet", contains="oligoSnpSet")
setClass("AnnotatedSnpCallSet", contains="SnpCallSet")
setClass("AnnotatedSnpCopyNumberSet", contains="SnpCopyNumberSet")

##no longer using oligo's definition of oligoSnpSet
##setClass("oligoSnpSet", contains="eSet")
setAs("AnnotatedSnpSet", "oligoSnpSet",
      function(from, to){
        warning("AnnotatedSnpSet DEPRECATED")
        new("oligoSnpSet",
            copyNumber=copyNumber(from),
            cnConfidence=cnConfidence(from),
            calls=calls(from),
            callsConfidence=callsConfidence(from),
            experimentData=experimentData(from),
            featureData=featureData(from),
            phenoData=phenoData(from),
            annotation=annotation(from))
      })

setAs("AnnotatedSnpCallSet", "SnpCallSet",
      function(from, to){
        warning("AnnotatedSnpCallSet DEPRECATED")
        new("SnpCallSet",
            calls=calls(from),
            callsConfidence=callsConfidence(from),
            experimentData=experimentData(from),
            featureData=featureData(from),
            phenoData=phenoData(from),
            annotation=annotation(from))
      })

setAs("AnnotatedSnpCopyNumberSet", "SnpCopyNumberSet",
      function(from, to){
        warning("AnnotatedSnpCopyNumberSet DEPRECATED")
        new("SnpCopyNumberSet",
            copyNumber=copyNumber(from),
            cnConfidence=cnConfidence(from),
            experimentData=experimentData(from),
            featureData=featureData(from),
            phenoData=phenoData(from),
            annotation=annotation(from))
      })


