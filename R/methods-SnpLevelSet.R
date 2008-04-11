##methods that I want to work for all the SNP Classes
setMethod("alleleA", "SnpLevelSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "SnpLevelSet", function(object) alleleB(featureData(object)))
setMethod("dbSnpId", "SnpLevelSet", function(object) dbSnpId(featureData(object)))

##setMethod("enzyme", "SnpLevelSet", function(object) enzyme(featureData(object)))
setMethod("fragmentLength", "SnpLevelSet", function(object) fragmentLength(fData(object)))

setMethod(".getX", "SnpLevelSet", function(object, ...) position(object))
setMethod(".getY", "SnpLevelSet", function(object, op, ...){
	if("copyNumber" %in% ls(assayData(object)) | "ratio" %in% ls(assayData(object))){
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
          



