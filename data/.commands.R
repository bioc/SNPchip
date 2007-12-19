##library(SNPchip) ##old version
##data(annSnpset)
##
##sample.snpset <- new("oligoSnpSet",
##                     calls=calls(annSnpset),
##                     callsConfidence=callsConfidence(annSnpset),
##                     copyNumber=copyNumber(annSnpset),
##                     cnConfidence=cnConfidence(annSnpset),
##                     annotation=annotation(annSnpset),
##                     featureData=featureData(annSnpset),
##                     phenoData=phenoData(annSnpset),
##                     experimentData=experimentData(annSnpset),
##                     annotation="pd.mapping50k.xba240")
##save(sample.snpset, file="~/projects/software/SNPchip/data/sample.snpset.RData", compress=TRUE)

##fD <- featureData(sample.snpset)
tmp <- new("AnnotatedDataFrame",
           data=data.frame(row.names=featureNames(sample.snpset)),
           varMetadata=data.frame(labelDescription=character()),
           dimLabels=c("featureNames", "featureColumns"))


tmp <- new("AnnotatedDataFrame")
tmp <- new("AnnotatedDataFrame",
           data=data.frame(),
           varMetadata=data.frame(),
           dimLabels=c("featureNames", "featureColumns"),
           
           
