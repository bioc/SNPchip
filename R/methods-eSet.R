setMethod("addFeatureData", "eSet", 
          function(object, path=NULL){
            ##This code will be removed as soon as annotation packages are available
            if(sum(annotation(object) == "mapping10k" |
                   annotation(object) == "mapping50kHind240" |
                   annotation(object) == "mapping50kXba240" |
                   annotation(object) == "mapping50kHind240:mapping50kXba240"|
                   annotation(object) == "mapping100k" |
                   annotation(object) == "mapping250kNsp" |
                   annotation(object) == "mapping250kSty" |
                   annotation(object) == "mapping250kNsp:mapping250kSty" |
                   annotation(object) == "mapping500k") < 1){
              stop("Annotation is only provided for the following Affymetrix platforms at this time:  mapping10k,  mapping50kHind240, mapping50kXba240, mapping50kHind240:mapping50kXba240, mapping100k (soon to be deprecated), mapping250kNsp, mapping250kSty, mapping250kNsp:mapping250kSty, mapping500k (soon to be deprecated)")
            }
            if(is.null(path)){
              print(paste("Loading annotation data from", path))
              annotation <- .getAnnotation(annotation(object))
            }
            if(!is.null(path)){
              if(annotation(object) == "mapping10k"){
                load(paste(path, annotation(object), ".rda", sep=""))
                annotation <- mapping10k
              }
              if(annotation(object) == "mapping50kHind240"){
                load(paste(path, "mapping50kHind240.rda", sep=""))
                annotation <-mapping50kHind240
              }
              if(annotation(object) == "mapping50kXba240"){
                load(paste(path, "mapping50kXba240.rda", sep=""))
                annotation <- mapping50kXba240
              }              
              if(annotation(object) == "mapping100k"){
                print("loading Hind...")
                load(paste(path, "mapping50kHind240.rda", sep=""))
                load(paste(path, "mapping50kXba240.rda", sep=""))
                annotation <- rbind(mapping50kHind240, mapping50kXba240)
              }
              if(annotation(object) == "mapping250kNsp"){
                load(paste(path, "mapping250kNsp.rda", sep=""))
                annotation <- mapping250kNsp
              }
              if(annotation(object) == "mapping250kSty"){
                load(paste(path, "mapping250kSty.rda", sep=""))
                annotation <- mapping250kSty
              }              
              if(annotation(object) == "mapping250kNsp:mapping250kSty"){
                load(paste(path, "mapping250kNsp.rda", sep=""))
                load(paste(path, "mapping250kSty.rda", sep=""))
                annotation <- rbind(mapping250kNsp, mapping250kSty)
              }
            }
            annotation <- as.matrix(annotation)
            annotation <- data.frame(annotation, stringsAsFactors=FALSE)
            annotation$Physical.Position <- as.integer(annotation$Physical.Position)            
            
            ##Exclude probes not on chromosomes 1-22,X, or Y
            chrom <- paste("chr", annotation$Chromosome, sep="")
            chrom[chrom == "chrX"] <- "chr23"
            chrom[chrom == "chrY"] <- "chr24"
            chr <- paste("chr", 1:24,sep="")
            annotation <- annotation[chrom %in% chr, ]

            ##Exclude probes not in the SnpSet object, and put the
            ##annotation object in the same order as the SnpSet object
            keep <- match(featureNames(object), annotation$Probe.Set.ID)
            ##NA's in keep means that the probe was not mapped to chr1 - 23, x, or Y
            ##Exclude these from the SnpSet
            object <- object[!is.na(keep), ]
            annotation <- annotation[keep[!is.na(keep)], ]
            if(!(identical(annotation$Probe.Set.ID, featureNames(object)))) stop("not matched")
            chrom <- chrom[keep[!is.na(keep)]]

            ###########################################################################
            ##Order the snps by chromosome and physical position
            chrom <- annotation$Chromosome
            chrom[chrom == "X"] <- "23"
            chrom <- as.numeric(chrom)
            ordering <- order(chrom, annotation$Physical.Position, decreasing=FALSE)
            ##keeping as factor really slows things down
            annotation$Chromosome <- paste("chr", annotation$Chromosome, sep="")
            object <- object[ordering, ]
            annotation <- annotation[ordering,]
            labeldescription <- colnames(annotation)
            vmd <- data.frame(labeldescription, row.names=colnames(annotation))
            colnames(vmd) <- "labelDescription"
            featureData(object) <- new("AnnotatedDataFrame",
                                       data=annotation, varMetadata=vmd)
            object
          })            
