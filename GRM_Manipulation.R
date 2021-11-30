### Script to manipulate GRMs:
# - Concatenate the GRM diagonals and off-diagonals
# - Filter out samples with diagonals above fixed thresholds


#Function to read a GRM
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}



# Concate all GRMs files in a directory into an object
GRM_alg0  <- list.files(pattern="*.grm.bin", path="./")
GRM_alg0  <- unlist(sapply(GRM_alg0, function(x) gsub(".grm.bin","",x)))
GRM_alg0 <- lapply(GRM_alg0,  function(x) ReadGRMBin(paste0("./", x)))

#Get the diagonals in a data frame, along with the bin number
alldiag <- data.frame(indiv = c(), diag = c(), bin = c())
for (i in c(1:8)) {
  test <- data.frame(indiv=as.vector(GRM_alg0[[i]]$id$V1), diag=GRM_alg0[[i]]$diag, bin = i)
  alldiag <- rbind(alldiag, test)
}

#Filter out samples that have a diagonal value larger than fixed thresholds (here outside of 0.7 and 1.3)
alldiagclean <- alldiag %>% filter(!diag > 1.3) %>% filter(!diag < 0.7)
alldiagclean %>% group_by(bin) %>% summarize(min = min(diag), max = max(diag))
length(which(table(alldiagclean$indiv) == 8))
indiv <- as.data.frame(c(which(table(alldiagclean$indiv) == 8)))
indiv <- as.data.frame(indiv)
indivlist <- data.frame(ID = rownames(indiv))

##### Extra code to manipulate a GRM elements directly (ie: set all diagonals to 1 or offdiagonals to 0) #####
size <- 4
for (i in c(1:2)){
    for (j in c(1:2)){
        GRM_read <- ReadGRMBin(paste0("./0", i, "_", j, "_GRM")) #Read GRM
        n <- length(GRM_read$diag) #Get sample size
        collapsed.grm <- vector(mode="numeric", length = n*(n+1)/2) #Create an empty vector of GRM values
        
        sum_i <- function(i){
            return(sum(1:i))
        }

        k <- sapply(1:n, sum_i) #Get index of diagonal elements
        collapsed.grm[k] <- GRM_read$diag #Assign the diagonals elements to their corresponding index
        collapsed.grm[k] <- 1 #Assign a value of 1 to all diagonal elements instead 
        collapsed.grm[-k] <- GRM_read$off #Assign the off-diagonal elements to their corresponding index
        collapsed.grm[-k] <- 0 #Assign a value of 0 to all off-diagonal elements instead 
        BinFileName <- as.character(paste0("./GRM_", i, "_", j, "_offdiagto0_diagonalsto1.grm.bin")) #Name the GRM
        BinFile <- file(BinFileName, "wb") #Create the file
        writeBin(con = BinFile, collapsed.grm, size = size) #Write the file, size = 4
        close(BinFile)
    }
}
