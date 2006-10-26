setMethod("initialize", "SnpCallSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   featureData = new("AnnotatedDataFrame"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix")) {
            callNextMethod(.Object,
                           assayData = assayDataNew(
                             calls = calls,
                             callsConfidence = callsConfidence),
                           featureData = featureData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })
setMethod("alleleA", "SnpCallSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "SnpCallSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "SnpCallSet", function(object) chromosome(featureData(object)))
setMethod("dbSnpId", "SnpCallSet", function(object) dbSnpId(featureData(object)))
setMethod("enzyme", "SnpCallSet", function(object) enzyme(featureData(object)))
setMethod("position", "SnpCallSet", function(object) position(featureData(object)))
setMethod("probeSetId", "SnpCallSet", function(object) probeSetId(featureData(object)))

setMethod("summary", signature(object = "SnpCallSet"),
          function(object, digits = 3, noCalls = FALSE, ...){
            het <- colMeans(ifelse(calls(object) == 2, 1, 0))
            hom <- colMeans(ifelse(calls(object) == 1 | calls(object) == 3, 1, 0))
            mat <- rbind(het, hom)
            mat
          })
