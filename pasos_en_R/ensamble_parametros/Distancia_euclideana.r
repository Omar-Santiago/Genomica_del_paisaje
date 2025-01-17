Distancia_euclideana <- function(liSNPs, param){ 
  liSNPs = liSNPs
  ##### Variables:
  # liSNPs = genlight object (package adegenet).
  # param = character string to label the stacks parameters or data used 
   
  require(adegenet)
  # 1) Get recovered sample-replicate pairs
  # get samples names
  samples<-indNames(liSNPs)
  
  # get succesfull replicate-sample pairs
  reps <- grep("r_",samples,value=TRUE) # get the replicates (ending with _r or _ir)
  samps <- match(sub("r_","_",reps),samples) # match against its sample (ie names w/o _r or _ir)
  samps<- samples[samps] # get names
  pairs<-cbind(reps,samps) # put them side by side 
  pairs<- pairs[rowSums(is.na(pairs)) < 1,] # remove rows that cointain NA to get only succesful replicate-sample pairs
  npairs <-nrow(data.frame(pairs))
  
  # 2) Estimate Euclidean distances between sample pairs 
  y <- numeric(0)
  for(i in 1:npairs) {
  
  srpair<-match(pairs[i,], samples) # 
  srpair <- liSNPs[srpair,] 
  indNames(srpair) # sample - replicate pair, first part of the name should be equal
  
  matriz<-dist(srpair, method = "euclidean")
  Eucl.dist<-matriz[1]

  # And put results toghether by sample-replicate pair
  pair <- paste(indNames(srpair)[2],"-",indNames(srpair)[1], sep="") # to generate a name for the pair
  Eclidean.distance<-cbind(param, pair, Eucl.dist)
  
  y<-rbind(y, Eclidean.distance)}

return(y)  
}
