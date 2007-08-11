##methods that I want to work for all the SNP Classes
setMethod("alleleA", "eSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "eSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "eSet", function(object) as.character(chromosome(featureData(object))))
setMethod("dbSnpId", "eSet", function(object) dbSnpId(featureData(object)))

setMethod("enzyme", "eSet", function(object) enzyme(featureData(object)))
setMethod("fragmentLength", "eSet", function(object) fragmentLength(fData(object)))

setMethod("getPar", "eSet", function(object, add.cytoband, ...){
  op <- switch(class(object),
               oligoSnpSet=new("ParSnpSet", ...),
               SnpCallSet=new("ParSnpCallSet", ...),
               SnpCopyNumberSet=new("ParSnpCopyNumberSet", ...))
  object <- object[!is.na(chromosome(object)), ]

  ##layout
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  N <- length(unique(chromosome(object)))
  S <- ncol(object)
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

  if(length(unique(chromosome(object))) > 10){
    op$alternate.xaxis.side <- TRUE
    side <- c(1, 3)[rep(1:2, (length(unique(chromosome(object)))/2) + 1)]
    side <- side[1:length(unique(chromosome(object)))]
    options(warn=-1)
    names(side) <- unique(chromosome(object))[order(as.numeric(unique(chromosome(object))))]
    op$xaxis.side <- side
  }
  m <- matrix(1:(S*N), nc=N, byrow=FALSE)
  w <- chromosomeAnnotation[unique(chromosome(object)), "chromosomeSize"]
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
  op
  ##set up defaults according to number of samples, chromosomes, position, etc.
})

setMethod(".getX", "eSet", function(object, ...) position(object))
          
setMethod("getSnpAnnotation", "eSet",
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

##Courtesy of Jason Ting
##setMethod("plotCytoband", "eSet",
##          function(object,
##                   cytoBand=NULL,
##                   xlim=NULL,
##                   cex.axis=0.8,
##                   xaxs="r",
##                   chromosome=NULL,
##                   main="",
##                   outer=TRUE,
##                   cytobandAxis=TRUE, ...){
##            
##            if(is.null(cytoBand)) {data(cytoband); cytoBand <- cytoband}
##            if(is.null(chromosome))  chrom <- unique(chromosome(object))  else chrom <- as.character(chromosome)[1]
##            cytoBand_chr <- cytoBand[as.character(cytoBand$chrom) == chrom,]
##            cytoBand_chr_p <- cytoBand_chr[grep("^p",as.character(cytoBand_chr$name)),]
##            cytoBand_chr_q <- cytoBand_chr[grep("^q",as.character(cytoBand_chr$name)),]
##            p.bands <- length(cytoBand_chr_p$chromEnd)
##            cut.left  <- c()
##            cut.right <- c()
##            ##  1st  band of arm or 1st  band after  "stalk"
##            ##  last band of arm or last band before "stalk"
##            for(band in 1:length(cytoBand_chr$chromEnd)) {
##              if (band==1)                             { cut.left[band] <- TRUE; cut.right[band] <- FALSE} else
##              if (band==p.bands)                       { cut.left[band] <- FALSE; cut.right[band] <- TRUE} else
##              if (band==(p.bands+1))                   { cut.left[band] <- TRUE; cut.right[band] <- FALSE} else
##              if (band==length(cytoBand_chr$chromEnd)) { cut.left[band] <- FALSE; cut.right[band] <- TRUE} else{
##                cut.left[band] <- FALSE; cut.right[band] <- FALSE
##              }
##            }
##            for(band in 1:length(cytoBand_chr$chromEnd)) {
##              if (as.character(cytoBand_chr$gieStain[band])=="stalk") {
##                cut.right[band-1] <- TRUE
##                cut.left[band] <- NA
##                cut.right[band] <- NA
##                cut.left[band+1] <- TRUE
##              }
##            }
##            ##When plotting subregions of a chromosome, this prevents the cytobands from extending beyond the subsetted object
##            ##exclude cytobands that end before the minimum plotting limits
##            include <- cytoBand_chr[, "chromEnd"] > xlim[1] & cytoBand_chr[, "chromStart"] < xlim[2]            
##            cytoBand_chr <- cytoBand_chr[include, ]
##            cut.left <- cut.left[include]
##            cut.right <- cut.right[include]            
##            plot(c(0, cytoBand_chr$chromEnd[length(cytoBand_chr$chromEnd)]),
##                 c(0, 2),
##                 xlim=xlim,
##                 type="n",
##                 xlab="",
##                 ylab="",
##                 axes=FALSE,
##                 xaxs=xaxs,
##                 main=main)
##            for (i in 1:length(cytoBand_chr$chromEnd)) {
##              start <- cytoBand_chr$chromStart[i]
##              end   <- cytoBand_chr$chromEnd[i]
##              delta = (end-start)/4
##              getStain <- function(stain){
##                switch(stain,
##                       gneg="grey100",
##                       gpos25="grey90",
##                       gpos50="grey70",
##                       gpos75="grey40",
##                       gpos100="grey0",
##                       gvar="grey100",
##                       stalk="brown3",
##                       acen="brown4")
##              }
##              color <- getStain(cytoBand_chr[i, "gieStain"])
####              stain <- as.character(cytoBand_chr$gieStain[i])
####              if (stain=="gneg") {
####                color <- "grey100"
####              } else if (stain=="gpos25") {
####                color <- "grey90"
####              } else if (stain=="gpos50") {
####                color <- "grey70"
####              } else if (stain=="gpos75") {
####                color <- "grey40"
####              } else if (stain=="gpos100") {
####                color <- "grey0"
####              } else if (stain=="gvar") {
####                color <- "grey100"
####              } else if (stain=="acen") {
####                color <- "brown4"
####              } else if (stain=="stalk") {
####                color <- "brown3"
####              } else {
####                color <- "white"
####              }
##              if (is.na(cut.left[i]) & is.na(cut.right[i])) {
##                ## this is a "stalk", do not draw box. Draw two vertival lines instead
##                delta <- (end-start)/3
##                lines(c(start+delta, start+delta), c(0,2), col=color)
##                lines(c(end-delta, end-delta), c(0,2), col=color)
##              } else if (cut.left[i] & cut.right[i]) {      # cut both ends
##                polygon(c(start, start+delta, end-delta, end, end, end-delta, start+delta, start),
##                        c(0.3, 0, 0, 0.3, 1.7, 2, 2, 1.7), col=color)
##              } else if (cut.left[i]) {              # cut left end only
##                polygon(c(start, start+delta, end, end, start+delta, start),
##                        c(0.3, 0, 0, 2, 2, 1.7), col=color)
##              } else if (cut.right[i]) {             # cut right end only
##                polygon(c(start, end-delta, end, end, end-delta, start),
##                        c(0, 0, 0.3, 1.7, 2, 2),col=color)
##              } else {
##                polygon(c(start, end, end, start),
##                        c(0, 0, 2, 2), col=color)
##              }
##            }
##            my.x <- (cytoBand_chr$chromStart+cytoBand_chr$chromEnd)/2
##            if(cytobandAxis){
##              axis(1, at=my.x, labels=as.character(cytoBand_chr$name), outer=outer, cex.axis=cex.axis,
##                   line=1, las=3)
##            }
##          })

##Call par only in the outer layer
##plots 1 sample and 1 chromosome
setMethod("plotSnp", "eSet",
          function(object, op, ...){
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
            
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
            options(warn=-1)
            objList <- objList[order(as.numeric(names(objList)))]
            options(warn=1)

            ##set graphical parameters for all plots
            par(allPlots(op))
            for(i in 1:length(objList)){
              if(i == 1) par(yaxt="s") else par(yaxt="n")
              .plotChromosome(objList[[i]], op=op)
            }
            mtext(op$ylab, side=op$side.ylab, outer=op$outer.ylab, las=3, cex=op$cex.ylab, line=op$line.ylab)
            mtext(op$xlab, side=op$side.xlab, outer=op$outer.xlab, cex=op$cex.xlab, line=op$line.xlab)            
          })


setMethod(".plotChromosome", "eSet",
          function(object, op, ...){
            cytoband <- .getCytoband(object, op)
            op$xlim <- .getXlim(object, op)
            for(j in 1:ncol(object)){
              .plot(object[, j], op=op)
              .drawCentromere(object[, j], op)              
              .drawCytobandWrapper(S=ncol(object), cytoband=cytoband, op=op, j=j)
              .drawXaxis(object=object, op=op, j=j)
            }
          })


setMethod("position", "eSet", function(object) position(featureData(object)))

