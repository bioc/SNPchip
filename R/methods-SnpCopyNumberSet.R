setMethod(".getY", "SnpCopyNumberSet", function(object, op, ...){
  y <- copyNumber(object)
  y[y < op$ylim[1]] <- op$ylim[1]
  y[y > op$ylim[2]] <- op$ylim[2]
  y
})

