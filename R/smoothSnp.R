#create object of class snpscan with smoothed copynumbers and smoothed loh calls
#LOH calls are smoothed by setting homozygous = 0 and heterozygous = 1
setMethod("smoothSnp", "AnnotatedSnpSet",
          function(chromosomes, object, samples,
                   span = 1/10,
                   method = "loess",
                   imputeNoCalls = TRUE,
                   verbose = TRUE){
            chromosomes <- paste("chr", chromosomes, sep = "")
            chromosomes[chromosomes == "chr23"] <- "chrX"
            chromosomes[chromosomes == "chr24"] <- "chrY"
            object <- object[chromosome(object) %in% chromosomes, samples]
            chrAnn <- chromosomeAnnotation(object)[chromosomes,]

            ##convert homozygous to 0 and heterozygous to 1
            calls(object) <- ifelse(calls(object) == 2, 1, 0)            
            smoothChromosome <- function(obj, span){
              
              loessX <- function(X, location, span){
                fit <- loess(X ~ location, span = span)$fitted
                return(fit)
              }

              cn.smooth <- apply(copyNumber(obj), 2, loessX, position(obj),
                                 span = span)
              call.smooth <- apply(calls(obj), 2, loessX,
                                   location = position(obj), span = span)

              copyNumber(obj) <- cn.smooth
              calls(obj) <- call.smooth
              obj
            }
            obj <- as(object, "AnnotatedSnpSetList")
            obj.list <- snpSetList(obj)
            obj@snpSetList <- lapply(obj.list, smoothChromosome, span = span)

            obj.new <- as(obj, "AnnotatedSnpSet")
            return(obj.new)
          })
