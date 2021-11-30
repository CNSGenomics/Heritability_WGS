#LD-partitioning

# This script was used to split each bim file (from different MAF bins) into ,4 list of variants grouped by LD score.
#It uses a *.score.ld file from GCTA as input. 

library(data.table)
library(tidyverse)
setwd("./{work_dir}")

#Load all BIM files
bimlist = list.files(pattern="./*.bim$")

#Read the LD file
allsnp <- fread("./{LD_file}.score.ld", h=T, select = c(1,8))


for (prefix in bimlist){
	lds_seg <- fread(prefix, select = 2, h=F)
	lds_seg_snp <- merge(lds_seg, allsnp, by.x="V2", by.y="SNP", all.x=F, all.y=F)
	lds_seg_snp <- lds_seg_snp %>%
		mutate(rank=percent_rank(ldscore_SNP))
	write.table(lds_seg_snp %>%
		 			filter(rank < 0.25) %>%
					select(V2),
				paste0(prefix,"_snp_group1.txt"), row.names=F, quote=F, col.names=F)
	write.table(lds_seg_snp %>%
					filter(rank >= 0.25 &  rank < 0.50) %>%
					select(V2),
				paste0(prefix,"_snp_group2.txt"), row.names=F, quote=F, col.names=F)
	write.table(lds_seg_snp %>%
					filter(rank >= 0.50 &  rank < 0.75) %>%
					select(V2),
				paste0(prefix,"_snp_group3.txt"), row.names=F, quote=F, col.names=F)
	write.table(lds_seg_snp %>%
					filter(rank >= 0.75) %>%
					select(V2),
				paste0(prefix,"_snp_group4.txt"), row.names=F, quote=F, col.names=F)
	
}