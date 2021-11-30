#### This script was used to perform QC on phenotypes (height and BMI) from the files provided with each TOPMed cohort.
# Each phenotype file across each cohort was first merged together in a single large phenotype file, {pheo_file}.phen.

# The structure of this {pheo_file}.phen file is 
# [1] "unique_subject_key"       "sample.id"                "subject_id"               "consent"                  "topmed_phs"               "study"                   
# [7] "topmed_project"           "topmed_phase"             "funding"                  "CCDG_ID"                  "seq_center"               "geno.cntl"               
#[13] "TRIO.dups"                "sex"                      "Xchr.anom"                "Ychr.anom"                "LDR"                      "exclude"                 
#[19] "SUBJECT_ID"               "topmed_study"             "height_baseline_1"        "age_at_height_baseline_1" "weight_baseline_1"        "age_at_weight_baseline_1"
#[25] "bmi_baseline_1"           "age_at_bmi_baseline_1" 

#Load the cleaned phenotypes available with genotype available (already merged by unique_key_id)
allpheno <- fread("{pheo_file}.phen", h=T)

#Load the list of unrelated europeans from PCA analysis and GRM --grm-singleton. THe output is a list of sample IDs
EURu <- fread("{Unrelated_EUR_file}.singleton.txt", h=F, select = c(1))

#Keep only the unrelated europeans in the phenotype file, remove unnecessary columns and samples aged < 18
unrel_eur <- allpheno[which(allpheno$sample.id %in% EURu$V1),]
unrel_eur <- unrel_eur %>% 
  filter(exclude == "FALSE") %>%
  select(-c(consent, topmed_project, topmed_phase, CCDG_ID, geno.cntl, Xchr.anom, Ychr.anom, exclude)) %>%
  filter(age_at_height_baseline_1 > 18) %>%
  drop_na()

#Proceed to phenotype standardization : group by study and sex, regress on age, scale residuals to get a variance of 1, and also created RINT height and BMI traits.

scale_this <- function(x){
  (x / sd(x))
}

unrel_eur <- unrel_eur %>% 
group_by(sex, study) %>% 
  mutate(height_resid = residuals(lm(height_baseline_1 ~ age_at_height_baseline_1)),
         bmi_resid = residuals(lm(bmi_baseline_1 ~ age_at_bmi_baseline_1))) %>%
  mutate(height_resid = scale_this(height_resid),
         bmi_resid = scale_this(bmi_resid)) %>% 
  mutate(height_resid_RINT =  qnorm((rank(height_resid, na.last = "keep")-0.5)/sum(!is.na(height_resid))),
         bmi_resid_RINT =  qnorm((rank(bmi_resid, na.last = "keep")-0.5)/sum(!is.na(bmi_resid)))) %>%
  ungroup() 

#Save files

write.table(unrel_eur %>% mutate(IID = sample.id) %>% select(IID, sample.id, height_resid),
            "./Height_EURu.phen", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(unrel_eur %>% mutate(IID = sample.id) %>% select(IID, sample.id, bmi_resid),
            "./BMI_EURu.phen", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(unrel_eur %>% mutate(IID = sample.id) %>% select(IID, sample.id, height_resid_RINT),
            "./Height_RINT_EURu.phen", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(unrel_eur %>% mutate(IID = sample.id) %>% select(IID, sample.id, bmi_resid_RINT),
            "./BMI_RINT_EURu.phen", col.names = F, row.names = F, quote = F, sep = "\t")
