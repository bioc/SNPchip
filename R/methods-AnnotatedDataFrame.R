setMethod("alleleA", "AnnotatedDataFrame", function(object) pData(object)$allele_a)
setMethod("alleleB", "AnnotatedDataFrame", function(object) pData(object)$allele_b)
setMethod("chromosome", "AnnotatedDataFrame", function(object) pData(object)$chrom)
setMethod("dbSnpId", "AnnotatedDataFrame", function(object) pData(object)$dbsnp_rs_id)
setMethod("position", "AnnotatedDataFrame", function(object) pData(object)$physical_pos)

setMethod("position", "AnnotatedDataFrame", function(object) pData(object)$physical_pos)
setMethod("enzyme", "AnnotatedDataFrame", function(object) pData(object)$enzyme)
setMethod("fragmentLength", "AnnotatedDataFrame", function(object) pData(object)$fragment_length)

##Copied from Biobase
setMethod("selectSomeIndex",
          signature(object="data.frame"),
          function(object, maxToShow=5, byrow=TRUE, ...) {
              len <-
                if (byrow) dim(object)[[1]]
                else dim(object)[[2]]
              if (maxToShow < 3) maxToShow <- 3
              if (len > maxToShow) {
                  maxToShow <- maxToShow - 1
                  bot <- ceiling(maxToShow/2)
                  top <- len-(maxToShow-bot-1)
                  list(1:bot, "...", top:len)
              } else if (len >= 1) list(1:len, NULL, NULL)
              else list(NULL, NULL, NULL)
          })
