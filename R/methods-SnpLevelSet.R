##methods that I want to work for all the SNP Classes
setMethod("alleleA", "SnpLevelSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "SnpLevelSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "SnpLevelSet", function(object) as.character(chromosome(featureData(object))))
setMethod("dbSnpId", "SnpLevelSet", function(object) dbSnpId(featureData(object)))

setMethod("enzyme", "SnpLevelSet", function(object) enzyme(featureData(object)))
setMethod("fragmentLength", "SnpLevelSet", function(object) fragmentLength(fData(object)))

setMethod("getPar", "SnpLevelSet", function(object, add.cytoband, ...){
  op <- switch(class(object),
               oligoSnpSet=new("ParSnpSet", ...),
               SnpCallSet=new("ParSnpCallSet", ...),
               SnpCopyNumberSet=new("ParSnpCopyNumberSet", ...))
  object <- object[!is.na(chromosome(object)), ]

  ##layout
  chromosomeNames <- unique(chromosome(object))
  chromosomeNames[chromosomeNames == "X"] <- 23
  chromosomeNames[chromosomeNames == "Y"] <- 24
  chromosomeNames <- chromosomeNames[order(as.numeric(chromosomeNames))]
  chromosomeNames[chromosomeNames == "23"] <- "X"
  chromosomeNames[chromosomeNames == "24"] <- "Y"  
  N <- length(chromosomeNames)
  S <- ncol(object)  
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  op$heights <- rep(1, ncol(object))
  if(!missing(add.cytoband)) op$add.cytoband <- add.cytoband
  if(op$add.cytoband){
    if(!op$outer.cytoband){
      S <- S+1
      if(S > 3){
        op$heights <- c(op$heights, 1/(ncol(object)*1.5))
      }
      if(S == 2 | S == 3){
        op$heights <- c(op$heights, 1/10)
      }
    } 
  }
  if(N > 10){
    op$alternate.xaxis.side <- TRUE
    side <- c(1, 3)[rep(1:2, N/2 + 1)]
    side <- side[1:N]
    options(warn=-1)
    names(side) <- chromosomeNames
    op$xaxis.side <- side
  }
  m <- matrix(1:(S*N), nc=N, byrow=FALSE)
  w <- chromosomeAnnotation[chromosomeNames, "chromosomeSize"]
  op$widths <- w/min(w)
  op$mat <- m

  ######################################################################
  ##graphical parameters
  if("copyNumber" %in% ls(assayData(object))){
    if(min(copyNumber(object), na.rm=TRUE) > 0) op$log <- "y"
    ##by default, we use the same ylimit on all the plots
    
    ##could make plot specific by adding an option one.ylim (or
    ##something to that effect) and calculating ylim in .plot()
    op$ylim <- .calculateYlim(object, op)
  }
  def.op <- options(warn=-1)
  op$firstChromosome <- chromosomeNames[1]
  options(def.op)
  op
  ##set up defaults according to number of samples, chromosomes, position, etc.
})

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

##Call par only in the outer layer
##plots 1 sample and 1 chromosome
setMethod("plotSnp", "eSet",
          function(object, op, on.exit=TRUE, ...){
            old.par <- par(no.readonly=TRUE)
            if(on.exit)  on.exit(par(old.par))
            
            object=object[!is.na(chromosome(object)), ]            
            if(missing(op)){
              op <- getPar(object, ...)
            }            
            if(op$useLayout){
              layout(mat=op$mat,
                     widths=op$widths,
                     heights=op$heights,
                     respect=op$respect)
            }
            objList <- split(object, chromosome(object))
            def.op <- options(warn=-1)

            ##this removes user control of the ordering.  Could add
            ##chromosomeOrder as an argument to ParESet
            names(objList)[names(objList) == "X"] <- "23"
            names(objList)[names(objList) == "Y"] <- "24"
            objList <- objList[order(as.numeric(names(objList)))]
            names(objList)[names(objList) == "23"] <- "X"
            names(objList)[names(objList) == "24"] <- "Y"            
            options(def.op)

            ##set graphical parameters for all plots
            par(allPlots(op))
            for(i in 1:length(objList)){
              if(i == 1) par(yaxt="s") else par(yaxt="n")
              .plotChromosome(objList[[i]], op=op)
            }
            if(op$outer.ylab) mtext(op$ylab, side=op$side.ylab, outer=TRUE, las=3, cex=op$cex.ylab, line=op$line.ylab)
            mtext(op$xlab, side=op$side.xlab, outer=op$outer.xlab, cex=op$cex.xlab, line=op$line.xlab)            
          })


setMethod(".plotChromosome", "eSet",
          function(object, op, ...){
            cytoband <- .getCytoband(object, op)
            op$xlim <- .getXlim(object, op)
            for(j in 1:ncol(object)){
              .plot(object[, j], op=op)
              .drawYaxis(object=object, op=op)
              .drawCentromere(object[, j], op)              
              .drawCytobandWrapper(S=ncol(object), cytoband=cytoband, op=op, j=j)
              .drawXaxis(object=object, op=op, j=j)
            }
          })


setMethod("position", "eSet", function(object) position(featureData(object)))

