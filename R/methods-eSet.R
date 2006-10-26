setMethod("addFeatureData", "eSet", 
          function(object, path=NULL){ #, snpAnnotation){

            if(sum(annotation(object) == "mapping10k" | annotation(object) == "mapping100k" |
                   annotation(object) == "mapping500k") < 1){
              stop("Annotation slot should be 'mapping10k', 'mapping100k', or 'mapping500k'")
            }
            if(is.null(path)){
              print(paste("Loading annotation data from", path))
              mapping <- .getAnnotation(annotation(object))
            }
            if(!is.null(path)){
              if(annotation(object) == "mapping10k"){
                load(paste(path, annotation(object), ".rda", sep=""))
                annotation <- mapping$annotation
              }
              if(annotation(object) == "mapping100k"){
                print("loading Hind...")
                load(paste(path, "mapping50kHind240.rda", sep=""))
                mappingHind <- mapping$annotation
                print("loading Xba...")
                load(paste(path, "mapping50kXba240.rda", sep=""))
                mappingXba <- mapping$annotation
                annotation <- rbind(mappingHind, mappingXba)
              }
              if(annotation(object) == "mapping500k"){
                load(paste(path, "mapping250kNsp.rda", sep=""))
                mappingNsp <- mapping500kNsp$annotation
                load(paste(path, "mapping250kSty.rda", sep=""))
                mappingSty <- mapping500kSty$annotation
                annotation <- rbind(mappingNsp, mappingSty)
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
