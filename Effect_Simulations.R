#### Script to add an effect to a phenotype simulated with GCTA.
# The input files are the .phen ones from GCTA. Example of a command line to simulate a phenotype:
$GCTA \
  --bfile ${bed_file} \
  --simu-qt \
  --simu-hsq 0.8 \
  --out ./${outfile} \
  --simu-causal-loci ${Causal_variants}.snplist \
  --threads "$ncpu"


### Add effect to a subset of samples then create a trait with a RINT
UKB_Rotherham <- data.frame(V1 = as.character(c(1:100))) #Example with 100 samples with an added effect among the whole set 

for (effect in c(0.5,1,2)) {  #Add a given effect to the phenotype
  df1 <- fread(paste0("${phenotype_file}.phen"), h=F)
  phen1 <- df1 %>%
             mutate(V3 = scale(V3)) %>%
             mutate(V3 = ifelse(V1 %in% UKB_Rotherham$V1, V3 + effect, V3)) %>%
             mutate(V3_RINT = qnorm((rank(V3, na.last = "keep")-0.5)/sum(!is.na(V3))))
  
  write.table(phen1 %>% select(V1, V2, V3), "${pheno_effect}.phen", col.names = F, row.names = F, quote=F, sep="\t")
  write.table(phen1 %>% select(V1, V2, V3_RINT), "${pheno_effect_RINT}.phen", col.names = F, row.names = F, quote=F, sep="\t")
}

add_effect_RINT()
