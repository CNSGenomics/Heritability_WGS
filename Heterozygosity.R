#### This script removes samples beyond 3 standard deviations based on the sample heterozygosity rate.
# It takes as imput a .het file (output from plink --het), concatenated across all bins (with the bin in column "Bin")


clean_het <- function(df){
  dfindiv <- df %>% group_by(Bin) %>% 
    filter(Rate > mean(Rate) + 3*sd(Rate) | Rate < mean(Rate) - 3*sd(Rate)) %>% ungroup() %>% distinct(IID) 
  df_clean <- df %>% filter(!IID  %in% dfindiv$IID)
  return(df_clean)
} #Function that removes 3 sd samples in each bin


list_het <- list.files(pattern="0*.het", path="./")
list_het_read <- lapply(list_het,  function(x) fread(x,h=T))
for (i in c(1:length(list_het))) {
  list_het_read[[i]]$Rate <- (list_het_read[[i]][,5] - list_het_read[[i]][,3]) / list_het_read[[i]][,5] #Compute the heterozygosity rate
  list_het_read[[i]]$Bin <- i
}
recap_het <- bind_rows(list_het_read)


#Here with 4 round of heterozygosity QC, can be adujsted.
het_14_bins_1 <- clean_het(het_14_bins)
het_14_bins_2 <- clean_het(het_14_bins_1)
het_14_bins_3 <- clean_het(het_14_bins_2)
het_14_bins_4 <- clean_het(het_14_bins_3)

write.table(het_14_bins_4 %>% select(IID) %>% distinct(), "./${list_distinct_sample_hetQCed}.indivlist", col.names = F, row.names = F, quote = F, sep = "\t")