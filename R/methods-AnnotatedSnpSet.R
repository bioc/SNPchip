###########################################################################
##updating deprecated classes
setMethod("updateObject", "AnnotatedSnpSet",
          function(object, ..., verbose=FALSE){
            new("oligoSnpSet",
                calls=calls(object),
                callsConfidence=callsConfidence(object),
                copyNumber=copyNumber(object),
                cnConfidence=cnConfidence(object),
                experimentData=experimentData(object),
                phenoData=phenoData(object),
                featureData=featureData(object))
          })
