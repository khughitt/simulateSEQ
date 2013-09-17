#-------------------- Rna-seq data ----------------------------#
# Generate a one sample based on:
# 1. natureGeneProfile
# 2. dispersion
# 3. s
geneProfile <- function(natureGeneProfile, dispersion, s){
  # set the gamma distribution parameters to match the biological expectation
  shape <- 1 / dispersion 
  rate <- shape / natureGeneProfile
  
  # draw gene counts from biolgical population
  n <- length(natureGeneProfile)
  lambda <- rgamma(n=n, shape=shape, rate=rate)
  
  # library prep. (ie. scale by size factor)
  libAdjCounts <- lambda * s
  
  # observed counts after "sequencing and alignment"
  obsGeneProfile <- rpois(n=n, lambda=libAdjCounts)  
  
  obsGeneProfile
}

# Generate n samples from a given:
# 1. natureGeneProfile
# 2. dispersion
# 3. libSizeFactors 
# where n is length(libSizeFactors)
generateBioSmpls <- function(natureGeneProfile, 
                             libSizeFactors, cond="smpl.", 
                             disp = 0.2,
                             trendDisp=FALSE){
  
  if(trendDisp){
    # trended biological coefficient of variation squared
    alphaf <- function(x) 0.2 + 0.8/x
    dispersion <- alphaf(natureGeneProfile) 
  }else{  
    dispersion <- disp
  }
  
  # generate a profile for each library size adjustment factor
  hold <- NULL
  i <- 0
  smpl <- paste(cond, 1:length(libSizeFactors), sep="")
  
  for(s in libSizeFactors){
    i <- i+1
    hold[[smpl[i]]] <- geneProfile(natureGeneProfile, dispersion, s)
  }
  
  res <- as.data.frame(hold)
  res
}

#-------------------- Microarray data -------------------------#
# Generate a one sample based on:
# 1. natureGeneProfile: mean
# 2. dispersion: sd
generateArray <- function(natureGeneProfile, dispersion){
   dat <- rnorm(length(natureGeneProfile), natureGeneProfile, dispersion)
   dat
}











