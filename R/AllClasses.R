setClass("HmmPredict", contains="SnpLevelSet",
	 representation(states="character",
			breakpoints="data.frame"))
setClassUnion("NULLorHmmPredict", c("NULL", "HmmPredict"))
setClass("ParESet",
	 representation(snpPar="list",
			snpset="SnpLevelSet",
			hmmPredict="NULLorHmmPredict"))

setClass("ParSnpCopyNumberSet", contains="ParESet")
setClass("ParSnpCallSet", contains="ParESet")
setClass("ParSnpSet", contains="ParSnpCopyNumberSet")

setValidity("ParESet", function(object){
	valid <- validObject(object@snpset)
	if(!valid) msg <- "snpset is not a valid object" else msg <- NULL
	valid <- valid && validObject(object@hmmPredict)
	if(!valid) msg <- c(msg, "hmmPredict is not a valid object") 
	return(msg)
})
setClass("RatioSnpSet", contains="SnpLevelSet")

###########################################################################
##DEPRECATED CLASSES
###########################################################################
##setClass("AnnotatedSnpSet", contains="oligoSnpSet")
##setClass("AnnotatedSnpCallSet", contains="SnpCallSet")
##setClass("AnnotatedSnpCopyNumberSet", contains="SnpCopyNumberSet")

###########################################################################
##Converting deprecated classes to current classes
###########################################################################
##setAs("AnnotatedSnpSet", "oligoSnpSet",
##      function(from, to){
##        warning("AnnotatedSnpSet DEPRECATED")
##        new("oligoSnpSet",
##            copyNumber=copyNumber(from),
##            cnConfidence=cnConfidence(from),
##            calls=calls(from),
##            callsConfidence=callsConfidence(from),
##            experimentData=experimentData(from),
##            featureData=featureData(from),
##            phenoData=phenoData(from),
##            annotation=annotation(from))
##      })
##setAs("AnnotatedSnpCallSet", "SnpCallSet",
##      function(from, to){
##        warning("AnnotatedSnpCallSet DEPRECATED")
##        new("SnpCallSet",
##            calls=calls(from),
##            callsConfidence=callsConfidence(from),
##            experimentData=experimentData(from),
##            featureData=featureData(from),
##            phenoData=phenoData(from),
##            annotation=annotation(from))
##      })
##setAs("AnnotatedSnpCopyNumberSet", "SnpCopyNumberSet",
##      function(from, to){
##        warning("AnnotatedSnpCopyNumberSet DEPRECATED")
##        new("SnpCopyNumberSet",
##            copyNumber=copyNumber(from),
##            cnConfidence=cnConfidence(from),
##            experimentData=experimentData(from),
##            featureData=featureData(from),
##            phenoData=phenoData(from),
##            annotation=annotation(from))
##      })


