##methods that I want to work for all the SNP Classes
setMethod("alleleA", "SnpLevelSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "SnpLevelSet", function(object) alleleB(featureData(object)))


##setMethod("chromosome", "SnpLevelSet", function(object) as.character(chromosome(featureData(object))))
setMethod("dbSnpId", "SnpLevelSet", function(object) dbSnpId(featureData(object)))

setMethod("enzyme", "SnpLevelSet", function(object) enzyme(featureData(object)))
setMethod("fragmentLength", "SnpLevelSet", function(object) fragmentLength(fData(object)))

##setMethod("getPar", "SnpLevelSet", function(object, add.cytoband, ...){


setMethod(".getX", "SnpLevelSet", function(object, ...) position(object))

setMethod(".getY", "SnpLevelSet", function(object, op, ...){
  if("copyNumber" %in% ls(assayData(object))){
    y <- copyNumber(object)
    y[y < op$ylim[1]] <- op$ylim[1]
    y[y > op$ylim[2]] <- op$ylim[2]
  } else {
    y <- calls(object)
    ##homozygous are 0's
    y[y == 3 | y == 1] <- 0
    ##heterozygous are 1's
    y[y == 2] <- 1
    y <- jitter(y, amount=0.05)
  }
  y
})
          
setMethod("getSnpAnnotation", "SnpLevelSet",
          function(object){
            warning("This method will be deprecated in the next release.")
            if(sum(annotation(object) == "pd.mapping50k.hind240" |
                   annotation(object) == "pd.mapping50k.xba240" |
                   annotation(object) == "pd.mapping250k.nsp" |
                   annotation(object) == "pd.mapping250k.sty") < 1){ 
              stop("Annotation is only provided for the following Affymetrix platforms at this time:  pd.mapping50k.hind240, pd.mapping50k.xba240,  pd.mapping250k.nsp, pd.mapping250k.sty")
            }
            require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package is not available"))
            getPositionInfo <- function(x){
              conn <- db(get(annotation(x)))
              snps <- paste("(", paste(featureNames(x), collapse="','"), ")", sep="'")
              sql <- paste("SELECT man_fsetid, dbsnp_rs_id, chrom, physical_pos, strand, allele_a, allele_b, fragment_length FROM featureSet WHERE man_fsetid IN",
                           snps, "ORDER BY man_fsetid")              
              dbGetQuery(conn, sql)
            }
            fD <- getPositionInfo(object)
            rownames(fD) <- fD$man_fsetid
            fD <- fD[, 2:8]
            if(!(identical(rownames(fD), featureNames(object)))){
              fD <- fD[match(featureNames(object), rownames(fD)), ]
              if(!(identical(rownames(fD), featureNames(object)))) stop("probes not matched")
            }
            getEnzyme <- function(x){
              switch(annotation(x),
                     pd.mapping50k.hind240="Hind",
                     pd.mapping50k.xba240="Xba",
                     pd.mapping250k.sty="Sty",
                     pd.mapping250k.nsp="Nsp")
            }
            fD$enzyme <- rep(getEnzyme(object), nrow(object))
            vmd <- data.frame(labelDescription=colnames(fD), row.names=colnames(fD))
            featureData <- new("AnnotatedDataFrame", data=fD, varMetadata=vmd)
            if(!validObject(featureData)) print("Not a valid AnnotatedDataFrame")
            featureData
          })



