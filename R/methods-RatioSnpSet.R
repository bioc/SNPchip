setMethod("copyNumber", "RatioSnpSet", function(object) ratio(object))
setMethod("ratio", "SnpLevelSet", function(object) assayDataElement(object, "ratio"))
setReplaceMethod("ratio", signature(object="SnpLevelSet", value="matrix"),
                 function(object, value){
			 assayDataElementReplace(object, "ratio", value)
		 })
