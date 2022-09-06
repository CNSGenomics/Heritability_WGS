########################################
# Script used to run GCTA/plink/bcftools for the analyses in the GREML-LDMS WGS 2022 paper
# It contains the most important commands for performing GREML-LDMS from WGS
# The Online Methods describe with greater details all the QC steps performed
# Additional R scripts are also available in this git
########################################

#### Variables
ncpu=12 
sample_list="/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/phenotypes/NTproBNP_14k.tsv" #List of samples with phenotypic information available

###### Create BED files and perform first filtering on SVM variants and other QC metrics


for i in  {1..22} ; do
    	bcftools view -m2 -M2  -Ou  -i'FILTER="PASS"' --threads "$ncpu" -S "$sample_list" --force-samples /data/project/Arora_lab/akhil/TOPMED/COMPLETE/GENOTYPES/freeze.10a.chr"$i".pass_only.gtonly.minDP10.vcf.gz | 
	bcftools  annotate --threads "$ncpu" -Ob -I +'%CHROM:%POS:%REF:%ALT' > /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/freeze10.14k.chr"$i".pass.bcf ;
	plink --bcf /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/freeze10.14k.chr"$i".pass.bcf --maf 0.0001 --allow-extra-chr --keep-allele-order --geno 0.05 --hwe 0.000001 --threads ${ncpu} -make-bed --out freeze10.14k.chr"$i".0.0001 ;
done


###### Merge BEDs ######
#${list_beds} contains the list of the autosomes (excepted chr 1) for merging

plink \
	--bfile ${BED_file_merged} \
	--merge-list ${list_beds} \
	--make-bed \
	--maf 0.0001 \
	--geno 0.05  \
	--hwe 0.000001 \
	--mind 0.05 \
	--out ${BED_file_merged_QC} \
	--threads ${ncpu}


###### PRS ######
#Simple PRS as a cohort QC check
#${scoreSNP} contains the effect sizes of the SNPs selected to construct the PRS

plink \
	--bfile ${BED_file_merged} \
	--threads ${ncpu}  \
	--score ${scoreSNP} \
	--out P


###### Computing IBD ######
#Computing IBD segments. For computational reasons, it is recomended to prune/select common SNPs for faster runtime.

king -b ${BED_file_merged_QC}.bed  \
	--ibdseg \
	--prefix ${BED_file_merged_QC_out_name} \
	--cpus ${ncpu} \
	--seglength 3


###### GRM ######
#Construct a GRM, extracting a list of variants (see 01_LD_bins.R to create LD bins from a bed of a MAF range)
#GRM constructed by part (one part computed per script on an array job)

i={1..99}
GCTA \
	--bfile ${BED_file_merged_QC} \
	--extract ${list_variants_LD_bin} \
	--make-grm-part 99 "$i" \
	--thread-num ${ncpu} \
	--out ${GRM_out} \
	--make-grm-alg 1


#Merge all GRM parts together

cat ${GRM_out}.part_99_*.grm.id > ${GRM_out}.grm.id
cat ${GRM_out}.part_99_*.grm.bin > ${GRM_out}.grm.bin
cat ${GRM_out}.part_99_*.grm.N.bin > ${GRM_out}.grm.N.bin


###### Extract unrelated samples ######
#Here we keep a relatedness threshold of 0.05

GCTA \
	--grm ${GRM_out} \
	--grm-cutoff 0.05 \
	--make-grm \
	--out ${GRM_out}_unrelated


###### Create a file containing multiple GRMs in a directory (need full path) ######

for i in *.grm.bin ; do readlink -f "$i"  | cut -d'.' -f1-2 >>  ${mgrm_file_path}; done


###### Filter variants wihin a MAF range ######

plink \
	--bfile ${BED_file_merged_QC} \
	--maf 0.0001 \
	--max-maf 0.001 \
	--make-bed \
	--threads ${ncpu} \
	--out ${BED_file_MAF}


###### SNP pruning ######
#Select independant SNPs based on a different thresholds from a set of variants
#Indep-pairwise parameters can be modified depending if the pruning is performed on common or rare variants 

i={1..22}
plink \
	--bfile ${BED_file_merged_QC} \
	--chr "$i" \
	--extract ${list_variants_bin} \
	--indep-pairwise 50 5 0.1 \
	--out ${out_indep_var}_chr"$i" \
	--threads ${ncpu}


###### Compute PCs ######
#Compute PC from GRM (GCTA) or from BED (plink2)

GCTA \
	--grm ${GRM_out} \
	--threads ${ncpu} \
	--pca 20 \
	--out ${PCA_out}


plink2 \
	--bfile ${BED_file_merged_QC} \
	--extract ${list_variants_bin} \
	--pca 20 approx \
	--out ${PCA_out} \
	--thread-num ${ncpu}


###### REML ######
#Run unconstrained REML from multiple GRMs (for GREML-LDMS), fitting PCs as quantitative covariates
#REML-no-lrt doesn't calculate the reduced model

GCTA \
	--reml \
	--mgrm ${mgrm_file_path} \
	--reml-no-constrain \
	--pheno ${phenotype_file} \
	--out ${REML_output_file} \
	--thread-num ${ncpu} \
	--qcovar ${PCA_out} \
	--reml-no-lrt





