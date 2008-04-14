setMethod("copyNumber", "RatioSnpSet", function(object) ratio(object))
setReplaceMethod("copyNumber", signature(object="RatioSnpSet", value="matrix"),
		 function(object, value){
			 ratio(object) <- value
			 object
		 })
setMethod("cnConfidence", "RatioSnpSet", function(object) ratioConfidence(object))
setMethod("ratioConfidence", "RatioSnpSet", function(object) assayDataElement(object, "ratioConfidence"))
setReplaceMethod("cnConfidence", signature(object="RatioSnpSet", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "ratioConfidence", value)
	})
setReplaceMethod("ratioConfidence", signature(object="RatioSnpSet", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "ratioConfidence", value)
		 })
setMethod("ratio", "SnpLevelSet", function(object) assayDataElement(object, "ratio"))
setReplaceMethod("ratio", signature(object="SnpLevelSet", value="matrix"),
                 function(object, value){
			 assayDataElementReplace(object, "ratio", value)
		 })
