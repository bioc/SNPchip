setMethod("calculateCopyNumber", signature(object = "SnpQSet"), function(object){
  A <- getA(object)
  cn <- rowMeans(A, dims = 2, na.rm = TRUE)
  cn
})
