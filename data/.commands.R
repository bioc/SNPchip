library(SNPchip) ##old version
data(annSnpset)

sample.snpset <- new("oligoSnpSet",
                 calls=calls(annSnpset),
                 callsConfidence=callsConfidence(annSnpset),
                 copyNumber=copyNumber(annSnpset),
                 cnConfidence=cnConfidence(annSnpset),
                 annotation=annotation(annSnpset),
                 featureData=featureData(annSnpset),
                 phenoData=phenoData(annSnpset),
                 experimentData=experimentData(annSnpset))
save(sample.snpset, file="~/projects/software/SNPchip/data/sample.snpset.RData", compress=TRUE)
