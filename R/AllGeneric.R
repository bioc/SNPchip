setGeneric("getSnpAnnotation",             function(object, ...) standardGeneric("getSnpAnnotation"))
setGeneric("alleleA", function(object) standardGeneric("alleleA"))
setGeneric("alleleB", function(object) standardGeneric("alleleB"))
setGeneric("chromosomeAnnotation", function(object) standardGeneric("chromosomeAnnotation"))
setGeneric("chromosomeAnnotation<-", function(object, value) standardGeneric("chromosomeAnnotation<-"))
setGeneric("dbSnpId", function(object) standardGeneric("dbSnpId"))
setGeneric("enzyme", function(object) standardGeneric("enzyme"))
setGeneric("fData", function(object) standardGeneric("fData"))
setGeneric("fData<-", function(object, value) standardGeneric("fData<-"))
setGeneric("fragmentLength", function(object) standardGeneric("fragmentLength"))
setGeneric("plotCytoband",        function(object, ...) standardGeneric("plotCytoband"))
setGeneric("smoothSnp", function(object, ...) standardGeneric("smoothSnp"))
setGeneric("plotSnp",        function(object,
                                      chromosomes,
                                      samples,
                                      snpId=c(NA, 1e6),
                                      oma=c(5, 4, 4, 0.5),
                                      mar=c(0.5, 0, 0.5, 0.2),
                                      width.right=NULL,
                                      showLayout=FALSE,
                                      plot=TRUE,
                                      col=c("royalblue", "red", "royalblue", "green3"),
                                      bg=rep("white", 4),
                                      cex=c(2, 3, 2, 2),
                                      pch=rep(".", 4),
                                      col.centromere=c("bisque", NA),
                                      bw=FALSE,
                                      lwd=2,
                                      col.axis="brown",
                                      cex.axis=1,
                                      cex.main=1,
                                      cex.lab=1,
                                      xlim=NULL,
                                      ylim=NULL,
                                      log="",
                                      xaxis.side=rep(c(1,3), length.out=length(chromosomes)),
                                      xaxs="i",                                      
                                      xaxt="n",
                                      yaxs="i",                                      
                                      yaxt="n",
                                      xlab="",
                                      ylab="",
                                      lab=c(2, 5, 7),
                                      adj=0,
                                      main="",
                                      legend=c(TRUE, TRUE),
                                      legend.panel=c(TRUE, FALSE),
                                      legend.location=c("topright", "bottomright"),
                                      legend.bty=c("n", "n"),
                                      legend.col=c("white", "white"),
                                      ncol=2,
                                      cex.legend=c(1,1),
                                      digits=3,
                                      pt.cex=cex.legend*1.5,
                                      bty="n",
                                      addCytoband=FALSE,
                                      height.cytoband=0.5, ...) standardGeneric("plotSnp"))


