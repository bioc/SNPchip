##setClass("SnpLayout", representation(mat="matrix",
##                                     widths="numeric",
##                                     heights="numeric",
##                                     respect="logical"))
setClass("ParESet", representation(snpPar="list"))
setClass("ParSnpCopyNumberSet", contains="ParESet")
setClass("ParSnpCallSet", contains="ParESet")
setClass("ParSnpSet", contains="ParSnpCopyNumberSet")

    


