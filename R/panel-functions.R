my.xypanel <- function(x, y,
		       x0, x1, chr.size,
		       col, border, coverage,
		       chr, show.coverage=TRUE,
		       max.y,
		       chromosomeAnnotation,
		       addCentromere=TRUE,
		       ..., subscripts){
	panel.grid(h=-1, v=10)
	panel.xyplot(x, y, ..., subscripts)
	h <- 0.75
	lrect(xleft=x0[subscripts],
	      xright=x1[subscripts],
	      ybottom=y-h/2,
	      ytop=y+h/2,
	      border=border[subscripts],
	      col=col[subscripts], ...)
	if(show.coverage)
		ltext(x, y,labels=coverage[subscripts], cex=0.6)
	##plot centromere
	if(addCentromere){
		chr <- unique(as.integer(as.character(chr)))
		coords <- chromosomeAnnotation[chr, 1:2]/1e6
		lrect(xleft=coords[1],
		      xright=coords[2],
		      ybottom=0,
		      ytop=max.y+h/2,
		      col="grey",
		      border="grey")
	}
}

xypanel <- function(x, y,
		    gt,
		    is.snp,
		    range,
		    col.hom="grey20",
		    fill.hom="lightblue",
		    col.het="grey20" ,
		    fill.het="salmon",
		    col.np="grey20",
		    fill.np="grey60",
		    show.state=TRUE,
		    cex.state=1,
		    col.state="blue",
		    cex.pch=0.3,
		    ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x[1], y[1], col="white", cex=cex.pch, ...) ## set it up, but don't plot
	is.snp <- is.snp[subscripts]
	if(!missing(gt)){
		gt <- gt[subscripts]
		hets.index <- which(gt == 2)
		hom.index <- which(gt == 1 | gt == 3)
		if(all(!c("col", "fill") %in% names(list(...)))){
			if(any(!is.snp))
				lpoints(x[!is.snp], y[!is.snp], col=col.np,
					fill=fill.np, cex=cex.pch, ...)
			if(length(hom.index) > 0)
				lpoints(x[hom.index], y[hom.index], col=col.hom,
					fill=fill.hom, cex=cex.pch, ...)

			if(length(hets.index) > 0)
				lpoints(x[hets.index], y[hets.index],
					col=col.het,
					fill=fill.het, cex=cex.pch, ...)
		}
	} else {
		lpoints(x[!is.snp], y[!is.snp], col=col.np,
			fill=fill.np, cex=cex.pch, ...)
		## use whatever col.hom to color SNPs
		lpoints(x[is.snp], y[is.snp], col=col.hom,
			fill=fill.hom, cex=cex.pch, ...)
	}
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
	if(show.state){
		## left justify the label to the start of the range
		y.max <- current.panel.limits()$ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1), cex=cex.state, col=col.state)
	}
}

xypanelBaf <- function(x, y,
		       gt,
		       baf,
		       is.snp,
		       range,
		       col.hom="grey20",
		       fill.hom="lightblue",
		       col.het="grey20" ,
		       fill.het="salmon",
		       col.np="grey20",
		       fill.np="grey60",
		       show.state=TRUE,
		       cex.state=1,
		       col.state="blue",
		       cex.pch=0.3,
		       ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x[1], y[1], col="white", cex=cex.pch, ...) ## set it up, but don't plot
	is.snp <- is.snp[subscripts]
	ylim <- current.panel.limits()$ylim
	y[y>ylim[2]] <- ylim[2]

	lpoints(x[!is.snp], y[!is.snp], col=col.np,
		fill=fill.np, cex=cex.pch, ...)
	## use whatever col.hom to color SNPs
	lpoints(x[is.snp], y[is.snp], col=col.hom,
		fill=fill.hom, cex=cex.pch, ...)
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
	if(show.state){
		## left justify the label to the start of the range
		y.max <- ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1), cex=cex.state, col=col.state)
	}
	b <- baf[subscripts]
	bnew <- rescale(b, ylim[1], ylim[1]+1.5)
	lpoints(x[is.snp], bnew[is.snp], cex=cex.pch, col="blue", ...)
}

prepanel.fxn <- function(x,y, chr.size, ..., subscripts){
	list(xlim=c(0, unique(chr.size[subscripts])), ylim=range(as.integer(as.factor(y[subscripts]))))
}