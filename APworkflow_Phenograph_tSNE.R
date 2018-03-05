
library(flowCore)
library(cytofkit)
library(flowMeans)
library(R.utils)
library(R.matlab)


# MAX TIME = 3 hours
maxtime_s = 9600	 

# Read iteration number
itmat=readMat("it.mat")
strit=paste("it",as.character(itmat),sep="")
print(strit)

datafr <- flowCore::read.FCS(paste("in_",strit,".fcs",sep=""), transformation = FALSE, truncate_max_range = FALSE)
data <- flowCore::exprs(datafr)

#------------
# Rphenograph
# (semiauto)
#------------
if (!file.exists(paste("Rpheno_",strit,".mat",sep=""))) {
  runtime1 <- system.time({
    tryCatch({
      Rpheno<-rep(NaN,dim(data)[1])
      Rpheno_aux <- withTimeout({
        Rphenograph(data,k=30)}, 
        timeout=maxtime_s, onTimeout = "error");
      Rpheno <- Rpheno_aux$membership
    },error= function(err){
      print(paste("MY_ERROR:  ",err))
      Rpheno<-rep(NaN,dim(data)[1])
    })
  })
  print(runtime1)
  #save(Rpheno, file = "Rpheno.RData")
  writeMat(paste("Rpheno_",strit,".mat",sep=""),Rpheno=Rpheno)
}



#------------
# tSNE
#------------
if (!file.exists(paste("tsne_xy_",strit,".mat",sep=""))) {
  runtime2 <- system.time({
    tryCatch({ 
      tsne_xy <- cytof_dimReduction(data, method = "tsne")
    },error= function(err){
      print(paste("MY_ERROR:  ",err))
    })
  })
  print(runtime2)
  #save(tsne_xy, file = "tsne_xy.RData")
  writeMat(paste("tsne_xy_",strit,".mat",sep=""),tsne_1=tsne_xy[,1],tsne_2=tsne_xy[,2])
}else{
  tsne_xy <- matrix(unlist(readMat(paste("tsne_xy_",strit,".mat",sep=""))), ncol = 2, byrow = FALSE)
}

