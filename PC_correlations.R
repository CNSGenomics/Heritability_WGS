#### This script can be used to estimate population stratification in a data set between 2 sets of PCs generated from the odds and the even chromosomes.
# - First, it generates correlations R^2 between sets of even and odd PC files (eigenvectors)
# - Second, it use the sample North and East birth coordinates in the UKB to create a distance matrix
# - Computes Moran's I index of autocorrelation
# - Compute the PC crossproduct and plot on the birth-coordinates map


library(tidyverse)
library(data.table)

#Computing PCs correlations

pc_even <- list.files(path="./Even/", pattern="eigenvec$")
pc_odd <- list.files(path="./Odd/", pattern="eigenvec$")

diag_recap <- data.frame(rsqpc=c(), bin=c())
N_pc <- 500

for (i in c(1:4)) {
  todd <- fread(paste0(pc_odd[i]), h=F, drop=c("V1", "V2"))
  teven <- fread(paste0(pc_even[i]), h=F, drop=c("V1", "V2"))
  todd <- todd[,1:N_pc]
  teven <- teven[,1:N_pc]
  diag_cor <- data.frame(rsqpc=apply(todd, 2, function(x){summary(lm(x ~ . , data=teven))$adj.r.squared}), bin=(i))
  diag_recap <-  rbind(diag_recap, diag_cor)
}

diag_recap_full <- diag_recap
diag_recap <- data.frame(rsqpc=c(), bin=c())

for (i in c(1:4)) {
  todd <- fread(paste0(pc_odd[i]), h=F, drop=c("V1", "V2"))
  teven <- fread(paste0(pc_even[i]), h=F, drop=c("V1", "V2"))
  todd <- todd[,1:N_pc]
  teven <- teven[,1:N_pc]
  diag_cor <- data.frame(rsqpc=apply(teven, 2, function(x){summary(lm(x ~ . , data=todd))$adj.r.squared}), bin=(i))
  diag_recap <-  rbind(diag_recap, diag_cor)
}

diag_recap_full <- diag_recap_full %>%
  mutate(PC=rep(c(1:N_pc), 2)) %>%
  mutate(meanrsq=(rsqpc + diag_recap$rsqpc)/2, variants = rep(c("Common", "Rare"),each=N_pc))
write.table(diag_recap_full, paste0("PCs_correlations_", N_pc), col.names=T, row.names=F, quote=F, sep="\t")


#Computing the samples distance matrix from east and north birth coordinates 
#File is IID, FID, #field

library(ape)

east <- fread("${east_coord}", h=F, select=c(1,3))
nord <- fread("${north_coord}", h=F, select=c(1,3))
coord <- merge(east, nord, by="V1")
setnames(coord, c("IID", "east_coord", "nord_coord"))
pc.dist.inv <- 1/as.matrix(dist(coord  %>% select(east_coord, nord_coord)))
diag(pc.dist.inv) <- 0
pc.dist.inv <- apply(pc.dist.inv, 2, function(x) ifelse(is.infinite(x), 0, x))

save(pc.dist.inv, file = "Distance_matrix_birth_coord_35k.RData")


#Compute Moran's I index for different MAF/LD bins and their associated Odd/Even eigenvector files

moran_res <- data.frame(maf = c(), ld=c(), pc=c(), moran=c())
moran_obs <-  function(maf, ld, pc){
  pc1 <- fread(paste0("./Even/0", maf, "_35k_Alg1_group", ld, "_Even.eigenvec"), h=F, select=c(3:152))
  pc2 <- fread(paste0("./Odd/0", maf, "_35k_Alg1_group", ld, "_Odd.eigenvec"), h=F, select=c(3:152))
  indiv <- fread(paste0("./Even/0", maf, "_35k_Alg1_group", ld, "_Even.eigenvec"), h=F, select=c(1))
  all <- as.data.frame(apply(pc1*pc2, 2, scale))
  all <-  as.data.frame(apply(all, 2, function(x){qnorm((rank(x, na.last = "keep")-0.5)/sum(!is.na(x)))}))

  #moran_res <- rbind(moran_res, data.frame(maf=maf, ld=ld, pc=c(1:20), moran=unlist(mclapply(allsub, function(x){ Moran.I(x, pc.dist.inv_sub)$observed}, mc.cores=6))))
  moran_res <- rbind(moran_res, data.frame(maf=maf, ld=ld, pc=c(1:150), moran=unlist(lapply(all, function(x){ Moran.I(x, pc.dist.inv)$observed}))))
  return(moran_res)
}


for(i in c(1:4)){
        for(j in c(1:2)){
                        moran_res <- moran_obs(i, j, k)
        }
}

write.table(moran_res, "Moran_results_full.OddEven", col.names=T, row.names=F, quote=F, sep="\t")

#Compute the PC crossproduct and plot on the coordinate maps

library(parallel)
library(ggplot2)
library(ape)

myIndex <- function(k,E,O,useAbs=FALSE){
  OE <- tcrossprod(O[k,],E[k,])
  if(useAbs){
    OE <- abs(OE)
  }
  return( 0.5* ( colMeans(OE) + rowMeans(OE) ) )
}
plot_popstrat2 <- function(maf, ld, pc){
  pc_col <- pc + 2
  maf_chr <- ifelse(maf > 3, "Common", "Rare")
  ld_chr <- ifelse(ld == 2, "HighLD", "LowLD")
  all <- cbind(coord, as.vector(indexMat[,k]))
  setnames(all, c("IID", "east_coord", "nord_coord", "pc_oe"))
  all <- all %>% mutate(pc_oe_scaled = scale(pc_oe)) %>% mutate(pc_oe_scaled = qnorm((rank(pc_oe_scaled, na.last = "keep")-0.5)/sum(!is.na(pc_oe))) )
  plot_ukb <- all  %>% filter(east_coord > 0.1 & nord_coord < 0.9) %>% ggplot(aes(x=east_coord, y=nord_coord, z=pc_oe_scaled)) +
    stat_summary_hex(fun = function(x) mean(x)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    labs(title=paste0("PC", pc,  "_Bin_", maf_chr, ld_chr, maf, "_", ld, " UKB WES"))
  ggsave(paste0(maf, "_", ld, "_PC", pc,  "_Odd_Even_Crossproduct_absFALSE", maf_chr, ld_chr, ".png"), plot_ukb, width = 8, height = 10, dpi = 300, units = "cm", scale=2.5)
}

for(i in c(1:4)){
  for(j in c(1:2)){
    east <- fread("UKB_Ex_East_Coord_35k", h=F, select=c(1,3))
    nord <- fread("UKB_Ex_North_Coord_35k", h=F, select=c(1,3))
    coord <- merge(east, nord, by="V1")
    coord <- coord %>% arrange(V1)
    pc_even_all <- fread(paste0("./Even/0", i, "_35k_Alg1_group", j, "_Even.eigenvec"), h=F)
    pc_odd_all <- fread(paste0("./Odd/0", i, "_35k_Alg1_group", j, "_Odd.eigenvec"), h=F)
    pc_even_all <- pc_even_all %>% arrange(V1) %>% select(-V1, -V2)
    pc_odd_all <- pc_odd_all %>% arrange(V1) %>% select(-V1, -V2)
    O <- as.matrix(pc_odd_all)
    E <- as.matrix(pc_even_all)
    N <- nrow(coord)
    indexMat <- do.call("rbind", mclapply(1:N, function(k) myIndex(k, E, O), mc.cores=12))
    for(k in c(1:5,50, 150)){
      plot_popstrat2(i, j, k)
    }
  }
}

