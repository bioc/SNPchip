#setMethod("summary", signature(object = "SnpCopyNumberSet"),
#          function(object, annotation, digits = 3, ...){
#            avg.cn <- colMeans(copyNumber(object))
#          })
#setMethod("copyNumber", "SnpCopyNumberSet", function(object) assayDataElement(object, "copyNumber"))
#setReplaceMethod("copyNumber", signature(object="SnpCopyNumberSet", value="matrix"),
#                 function(object, value) assayDataElementReplace(object, "copyNumber", value))
setMethod("initialize", "SnpCopyNumberSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   featureData = new("AnnotatedDataFrame"),
                   annotation = character(),
                   copyNumber = new("matrix"),
                   cnConfidence = new("matrix")) {
            callNextMethod(.Object,
                           assayData = assayDataNew(
                             copyNumber = copyNumber,
                             cnConfidence = cnConfidence),
                           featureData = featureData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })
setMethod("alleleA", "SnpCopyNumberSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "SnpCopyNumberSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "SnpCopyNumberSet", function(object) chromosome(featureData(object)))
setMethod("dbSnpId", "SnpCopyNumberSet", function(object) dbSnpId(featureData(object)))
setMethod("enzyme", "SnpCopyNumberSet", function(object) enzyme(featureData(object)))
setMethod("position", "SnpCopyNumberSet", function(object) position(featureData(object)))
##setMethod("probeSetId", "SnpCopyNumberSet", function(object) probeSetId(featureData(object)))
