setMethod("chromosome", "SnpCopyNumberSet",
          function(object){
            if(!("chromosome" %in% fvarLabels(object))){
              require("RSQLite") || stop("if 'chromosome' is not a varLabel in the featureData, the RSQLite package must be installed")
              fs <- featureNames(object)
              sql <- "SELECT man_fsetid, chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%'"              
              ##Check if two objects have been combined
              pkgs <- strsplit(annotation(object), ",")[[1]]
              if(length(pkgs) > 1){
                object2 <- object1 <- object
                annotation(object1) <- pkgs[1]
                annotation(object2) <- pkgs[2]

                tmp1 <- dbGetQuery(db(object1), sql)
                tmp2 <- dbGetQuery(db(object2), sql)
                tmp <- rbind(tmp1, tmp2)
              } else {
                tmp <- dbGetQuery(db(object), sql)
              }
              idx <- match(fs, tmp[["man_fsetid"]])
              chr <- tmp[idx, "chrom"]
            } else {
              chr <- featureData(object)$chromosome
            }
            return(chr)
          })

setMethod("position", "SnpCopyNumberSet",
          function(object){
            if(!("position" %in% fvarLabels(object))){
              require("RSQLite") || stop("If 'position' is not a varLabel in featureData, the RSQLite package must be installed")
              fs <- featureNames(object)
              sql <- "SELECT man_fsetid, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
              pkgs <- strsplit(annotation(object), ",")[[1]]
              if(length(pkgs) > 1){
                object2 <- object1 <- object
                annotation(object1) <- pkgs[1]
                annotation(object2) <- pkgs[2]                
              
                tmp1 <- dbGetQuery(db(object1), sql)
                tmp2 <- dbGetQuery(db(object2), sql)
                tmp <- rbind(tmp1, tmp2)
              } else{
                tmp <- dbGetQuery(db(object), sql)
              }
              idx <- match(fs, tmp[["man_fsetid"]])
              pos <- tmp[idx, "physical_pos"]
            } else {
              pos <- featureData(object)$position
            }
            return(pos)
          })



setMethod("getSnpAnnotation", "SnpCopyNumberSet",
          function(object){
            callNextMethod()
          })

