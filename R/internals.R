.getAnnotation <- function(chip){
  load10k <- function(){
    print("Fetching 10k annotation from http://biostat.jhsph.edu/...")
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping10k.rda")))
    return(mapping)
  }
  load100k <- function(){
#    print("Fetching 100k annotation from http://biostat.jhsph.edu/...")
    print("Fetching 50k Hind annotation from http://biostat.jhsph.edu/...")
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kHind240.rda")))
    mappingHind <- mapping$annotation
    print("Fetching 50k Xba annotation from http://biostat.jhsph.edu/...")    
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kXba240.rda")))
    mappingXba <- mapping$annotation
    mapping <- rbind(mappingHind, mappingXba)
    return(mapping)
  }              
  load500k <- function(){
    print("Fetching 100k annotation from http://biostat.jhsph.edu/...")                                
    try(load("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping250kNsp.rda"))
    mappingNsp <- mapping500kNsp$annotation
    try(load("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping250kSty.rda"))
    mappingSty <- mapping500kSty$annotation
    mapping <- rbind(mappingNsp, mappingSty)
    return(mapping)
  }
  switch(chip,
         mapping10k=load10k(),
         mapping100k=load100k(),
         mapping500k=load500k())
}




