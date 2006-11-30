.getAnnotation <- function(chip){
  if(chip == "mapping50kHind240:mapping50kXba240") chip <- "mapping100k"
  if(chip == "mapping250kNsp:mapping250kSty") chip <- "mapping500k"
  load10k <- function(){
    print("Fetching 10k annotation from http://biostat.jhsph.edu/...")
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping10k.rda")))
    return(mapping10k)
  }
  load50kXba <- function(){
    print("Retrieving 50k Xba annotation from http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/")    
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kXba240.rda")))
    mapping50kXba240
  }
  load50kHind <- function(){
    print("Retrieving 50k Hind annotation from http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/")
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kHind240.rda")))
    mapping50kHind240
  }
  load100k <- function(){
    print("Retrieving 50k Hind annotation from http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/")
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kHind240.rda")))
    print("Retrieving 50k Xba annotation from http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/")    
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping50kXba240.rda")))
    mapping100k <- rbind(mapping50kHind240, mapping50kXba240)
    return(mapping100k)
  }
  load250kNsp <- function(){
    print("Fetching 100k annotation from http://biostat.jhsph.edu/...")                                
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping250kNsp.rda")))
    mapping250kNsp
  }
  load250kSty <- function(){
    print("Fetching 100k annotation from http://biostat.jhsph.edu/...")                                
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping250kSty.rda")))
    mapping250kSty
  }  
  load500k <- function(){
    print("Fetching 100k annotation from http://biostat.jhsph.edu/...")                                
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping250kNsp.rda")))
    try(load(url("http://biostat.jhsph.edu/~iruczins/publications/sm/2006.scharpf.bioinfo/mapping/mapping250kSty.rda")))
    mapping500k <- rbind(mapping250kNsp, mapping250kSty)
    return(mapping500k)
  }
  switch(chip,
         mapping10k=load10k(),
         mapping50kXba240=load50kXba(),
         mapping50kHind240=load50kHind(),
         mapping100k=load100k(),
         mapping250kNsp=load250kNsp(),
         mapping250kSty=load250kSty(),
         mapping500k=load500k())
}




