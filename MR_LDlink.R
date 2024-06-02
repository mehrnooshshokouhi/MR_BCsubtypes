## The aim of this Rscript is to identify BCAC proxies for MR analyses ##

# Load required packages into R
library(TwoSampleMR)
library(openxlsx)
library(dplyr)
library(LDlinkR)

# Set working directory
setwd("/DATA/users/r.verdiesen/MR_BCAC/")


# load BCAC GWAS data into R 

lumA_data = TwoSampleMR::read_outcome_data(
  filename = "sumstat_BC/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",
  sep = " ",
  snp_col = "SNP.Onco",
  beta_col = "Luminal_A_log_or_meta",
  se_col = "Luminal_A_se_meta",
  eaf_col = "EAFcontrols.Onco", # use OncoArray info - because highest nr of participants
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  pval_col = "") # no p-value available; TwoSampleMR package infers p-values

head(lumA_data) # correct


# Make rsid variable
lumA_data$SNP_old = lumA_data$SNP

lumA_data$SNP = sapply(strsplit(lumA_data$SNP_old, 
                                split = ":"), "[", 1) # this does not result in an rsid for every entry;

head(lumA_data) # correct
lumA_data$SNP_old = NULL




## Define functions used later on for the selection of proxy SNPs ##

# 1. Function to make variables to indicate index SNPs and overlap with BCAC data
variables_proxies = function(proxies) {
  proxies %>%
    mutate(index_snp = case_when(
      RS_Number == query_snp ~ 1,
      RS_Number != query_snp ~ 0)) %>%
    mutate(in_BCAC = case_when(
      RS_Number %in% lumA_data$SNP ~ 1))
}


# 2. Function to select one BCAC proxy for each index SNP based on the highest R2
selecting_proxies_in_BCAC = function(proxies) {
  proxies %>%
    filter(! RS_Number %in% query_snp) %>% # only selects proxy SNPs that are not an index SNP themselves
    filter(in_BCAC == 1) %>% # only selects SNPs for which BCAC data is available
    filter(R2 >= 0.8) %>% # only selects SNPs that are a good proxy based on LD R2
    group_by(query_snp) %>% # perform filtering on R2 for each index SNP "group"
    filter(R2 == max(R2)) %>% # select SNPs with highest LD R2
    slice_head(n = 1) # only keep first entry in case there are multiple proxy SNPs with exactly the same LD R2 value
}


# 3. Function to select one BCAC proxy for each index SNP based on the highest R2, excluding indels;
# Only use this function if none of the index variants are indels!
selecting_proxies_in_BCAC_excl_indels = function(proxies) {
  proxies %>%
    filter(! RS_Number %in% query_snp) %>% # only selects proxy SNPs that are not an index SNP themselves
    filter(in_BCAC == 1) %>% # only selects SNPs for which BCAC data is available
    filter(R2 >= 0.8) %>% # only selects SNPs that are a good proxy based on LD R2
    filter(!grepl('-', Alleles)) %>% # filters out indels
    group_by(query_snp) %>% # perform filtering on R2 for each index SNP "group"
    filter(R2 == max(R2)) %>% # select SNPs with highest LD R2
    slice_head(n = 1) # only keep first entry in case there are multiple proxy SNPs with exactly the same LD R2 value
}



## Search proxies for each risk factor ##

# 1. Height

height_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_stats_height_sexcombined.txt.gz",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA_COJO",
  se_col = "SE_COJO",
  eaf_col = "Freq_Tested_Allele_in_HRS",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  pval_col = "P_COJO",
  samplesize_col = "N") # correct, 3290 SNPs

# Extract list with overlapping variants between two height and BCAC data frames
overlap_height_bcac = intersect(height_data$SNP,
                                lumA_data$SNP) # 2,884 SNPs in both data frames


# Extract SNPs that are in height data frame but NOT in BCAC data frame
proxy_needed_bcac_height = setdiff(height_data$SNP,
                                   lumA_data$SNP) # 406 height SNPs that do not overlap with BCAC

# Check if there are indels among the index variants
table(height_data$effect_allele.exposure);
table(height_data$other_allele.exposure) # no indels/triallelic snps als index variants; so use selecting_proxies_in_BCAC_excl_indels later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_height,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #


# Load list with height proxies into R  
proxies_height = read.table(file = "combined_query_snp_list.txt",
                            header = T,
                            sep = "\t",
                            row.names = NULL) # n = 605,308

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_height = setdiff(proxy_needed_bcac_height,
                              proxies_height$query_snp)
table(duplicated(not_searched_height)) # 13 SNPs; no duplicates

# Check data
head(proxies_height) # correct

table(duplicated(proxies_height$query_snp)) # correct; 393 + 13 = 406


# Perform proxy selection for height

# Add variables index_SNP and in_BCAC to data
proxies_height_2 = variables_proxies(proxies_height)

table(proxies_height_2$index_snp)
table(proxies_height_2$in_BCAC) # 295,386 proxies in BCAC


# Perform actual selection
proxies_height_final = selecting_proxies_in_BCAC_excl_indels(proxies_height_2) # n = 270 proxies
# Check if there is overlap between the index snps and selected proxies
check = intersect(proxies_height_final$query_snp,
                  proxies_height_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_height_final$RS_Number,
                   height_data$SNP) # this is also not the case!

rm(check, check2)


# Check if there are duplicates among the selected proxies
table(duplicated(proxies_height_final$RS_Number)) # there is one duplicated SNP id, which one?

print(proxies_height_final[which(duplicated(proxies_height_final$RS_Number)), c(1:3)]) # rs12423004
print(proxies_height_final[proxies_height_final$RS_Number == "rs12423004", c(1:3)]) # are there other proxies available for rs17032823?

print(proxies_height[proxies_height$query_snp=="rs17032823", c(1:3)]) # yes; lots!
head(proxies_height[proxies_height$query_snp=="rs17032823", c(1:9)]) 
# rs79544094 is a good candidate; R2 is also 1 and this SNP has not been selected as proxy for another query_snp

# First remove initial row for query_snp rs17032823
proxies_height_final_2 = proxies_height_final %>%
  filter(query_snp != "rs17032823")

# Then, add this the new proxy to proxies_height_final2
proxies_height_final_3 = rbind(proxies_height_final_2,
                               proxies_height[(proxies_height$query_snp == "rs17032823" & proxies_height$RS_Number == "rs79544094"),])

table(duplicated(proxies_height_final_3$RS_Number))

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_height_final_3$query_snp,
                  proxies_height_final_3$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_height_final_3$RS_Number,
                   height_data$SNP) # this is also not the case!

table(duplicated(proxies_height_final_3$RS_Number)) # correct; no duplicates left among selected proxies

rm(check, check2)


# Check R2 values
summary(proxies_height_final_3$R2) # correct, no R2 values below 0.8




# Save height proxy dataframe as .txt file
#  write.table(proxies_height_final_3,
#              file = "proxy_search/proxysnps_height_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # final proxy file saved at 20/12/2021


# Clean R environment; remove height related objects
rm(height_data,
   proxies_height,
   proxies_height_2,
   proxies_height_final,
   proxies_height_final_2,
   proxies_height_final_3,
   overlap_height_bcac,
   proxy_needed_bcac_height,
   not_searched_height)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 10:28h



# 2. Body mass index

bmi_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_stats_bmi_sexcombined.txt.gz",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA_COJO",
  se_col = "SE_COJO",
  eaf_col = "Freq_Tested_Allele_in_HRS",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  pval_col = "P_COJO",
  samplesize_col = "N") # correct, 941 SNPs


# Extract list with overlapping variants between two BMI and BCAC data frames
overlap_bmi_bcac = intersect(bmi_data$SNP,
                             lumA_data$SNP) # 835 SNPs in both data frames


# Extract SNPs that are in BMI data frame but NOT in BCAC data frame
proxy_needed_bcac_bmi = setdiff(bmi_data$SNP,
                                lumA_data$SNP) # 106 bmi SNPs that do not overlap with BCAC


# Check if there are indels among the index variants
table(bmi_data$effect_allele.exposure);
table(bmi_data$other_allele.exposure) # no indels snps als index variants; so use selecting_proxies_in_BCAC_excl_indels later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_bmi,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #


# Load list with bmi proxies into R  
proxies_bmi = read.table(file = "combined_query_snp_list.txt",
                         header = T,
                         sep = "\t",
                         row.names = NULL) # n = 170,659

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_bmi = setdiff(proxy_needed_bcac_bmi,
                           proxies_bmi$query_snp) 
table(duplicated(not_searched_bmi)) # 4 SNPs; no duplicates

# Check data
head(proxies_bmi) # correct

table(duplicated(proxies_bmi$query_snp)) # correct; 102 + 4 = 106


# Perform proxy selection for bmi

# Add variables index_SNP and in_BCAC to data
proxies_bmi_2 = variables_proxies(proxies_bmi)

table(proxies_bmi_2$index_snp)
table(proxies_bmi_2$in_BCAC) # 82,707 proxies available in BCAC


# Perform actual selection
proxies_bmi_final = selecting_proxies_in_BCAC_excl_indels(proxies_bmi_2) # n = 85 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_bmi_final$query_snp,
                  proxies_bmi_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_bmi_final$RS_Number,
                   bmi_data$SNP) # this is also not the case!

table(duplicated(proxies_bmi_final$RS_Number)) # no duplicates among selected proxies either

rm(check, check2)

# Check R2 values
summary(proxies_bmi_final$R2) # correct, no R2 values below 0.8


# Save BMI proxy dataframe as .txt file
#  write.table(proxies_bmi_final,
#              file = "proxy_search/proxysnps_bmi_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# Clean R environment; remove height related objects
rm(bmi_data,
   proxies_bmi,
   proxies_bmi_2,
   proxies_bmi_final,
   overlap_bmi_bcac,
   proxy_needed_bcac_bmi,
   not_searched_bmi)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 10:53h



# 2. Type 2 diabetes

t2d_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/index_t2d.txt",
  sep = "\t",
  snp_col = "rsid",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  pval_col = "P") # correct, 425 SNPs


# Extract list with overlapping variants between two T2D and BCAC data frames
overlap_t2d_bcac = intersect(t2d_data$SNP,
                             lumA_data$SNP) # 375 SNPs in both data frames


# Extract SNPs that are in T2D data frame but NOT in BCAC data frame
proxy_needed_bcac_t2d = setdiff(t2d_data$SNP,
                                lumA_data$SNP) # 50 t2d SNPs that do not overlap with BCAC

# Remove "-" from proxy_needed_bcac_t2d
proxy_needed_bcac_t2d_2 = proxy_needed_bcac_t2d[-30]

# Check if there are indels among the index variants
table(t2d_data$effect_allele.exposure);
table(t2d_data$other_allele.exposure) # no indels snps als index variants; so use selecting_proxies_in_BCAC_excl_indels later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_t2d_2,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #

# Load list with t2d proxies into R  
proxies_t2d = read.table(file = "combined_query_snp_list.txt",
                         header = T,
                         sep = "\t",
                         row.names = NULL) # n = 70,933

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_t2d = setdiff(proxy_needed_bcac_t2d,
                           proxies_t2d$query_snp) # 4 SNPs; no duplicates

# Check data
head(proxies_t2d) # correct

table(duplicated(proxies_t2d$query_snp)) # correct; 46 + 4 = 50


# Perform proxy selection for t2d

# Add variables index_SNP and in_BCAC to data
proxies_t2d_2 = variables_proxies(proxies_t2d)

table(proxies_t2d_2$index_snp)
table(proxies_t2d_2$in_BCAC) # 25,020 proxies in BCAC


# Perform actual selection
proxies_t2d_final = selecting_proxies_in_BCAC_excl_indels(proxies_t2d_2) # n = 29 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_t2d_final$query_snp,
                  proxies_t2d_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_t2d_final$RS_Number,
                   t2d_data$SNP) # this is also not the case!

table(duplicated(proxies_t2d_final$RS_Number)) # no duplicates among selected proxies either

rm(check, check2)


# Check R2 values
summary(proxies_t2d_final$R2) # correct, no R2 values below 0.8


# Save T2D proxy dataframe as .txt file
#  write.table(proxies_t2d_final,
#              file = "proxy_search/proxysnps_t2d_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# Clean R environment; remove height related objects
rm(t2d_data,
   proxies_t2d,
   proxies_t2d_2,
   proxies_t2d_final,
   overlap_t2d_bcac,
   proxy_needed_bcac_t2d,
   proxy_needed_bcac_t2d_2,
   not_searched_t2d)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 11:13h


# 4. Age at menarche

menarche_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/indexsnps_menarche.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA (y/allele)",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "Effect allele",
  other_allele_col = "Other allele",
  pval_col = "P",
  samplesize_col = "N") # correct, 389 SNPs

# Extract list with overlapping variants between age at menarche and BCAC data frames
overlap_menarche_bcac = intersect(menarche_data$SNP,
                                  lumA_data$SNP) # 345 SNPs in both data frames


# Extract SNPs that are in age at menarche data frame but NOT in BCAC data frame
proxy_needed_bcac_menarche = setdiff(menarche_data$SNP,
                                     lumA_data$SNP) # 44 SNPs that do not overlap with BCAC

# Check if there are indels among the index variants
table(menarche_data$effect_allele.exposure);
table(menarche_data$other_allele.exposure) # indels as index variants; so use selecting_proxies_in_BCAC later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_menarche,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #


# Load list with menarche proxies into R  
proxies_menarche = read.table(file = "combined_query_snp_list.txt",
                              header = T,
                              sep = "\t",
                              row.names = NULL) # n = 67,865

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_menarche = setdiff(proxy_needed_bcac_menarche,
                                proxies_menarche$query_snp) # 8 SNPs; no duplicates

# Check data
head(proxies_menarche) # correct

table(duplicated(proxies_menarche$query_snp)) # correct; 36 + 8 = 44


# Perform proxy selection for menarche

# Add variables index_SNP and in_BCAC to data
proxies_menarche_2 = variables_proxies(proxies_menarche)

table(proxies_menarche_2$index_snp) # 32 index SNPs
table(proxies_menarche_2$in_BCAC) # 41,302 proxies in BCAC

table(proxies_menarche_2$query_snp==proxies_menarche_2$RS_Number) # correct that there are 32 index snps instead of 36


# Perform actual selection
proxies_menarche_final = selecting_proxies_in_BCAC(proxies_menarche_2) # n = 24 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_menarche_final$query_snp,
                  proxies_menarche_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_menarche_final$RS_Number,
                   menarche_data$SNP) # this is also not the case!

table(duplicated(proxies_menarche_final$RS_Number)) # no duplicates among selected proxies either

rm(check, check2)


# Check R2 values
summary(proxies_menarche_final$R2) # correct, no R2 values below 0.8


# Save menarche proxy dataframe as .txt file
#  write.table(proxies_menarche_final,
#              file = "proxy_search/proxysnps_menarche_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# Clean R environment; remove height related objects
rm(menarche_data,
   proxies_menarche,
   proxies_menarche_2,
   proxies_menarche_final,
   overlap_menarche_bcac,
   proxy_needed_bcac_menarche,
   not_searched_menarche)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 11:32h


# 5. Age at menopause

menopause_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/new_menopause.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "Effect...19",
  se_col = "SE...20",
  eaf_col = "Effect Allele Freq.",
  effect_allele_col = "Effect Allele",
  other_allele_col = "Other Allele",
  pval_col = "P-value...21") # correct, 290 SNPs


# Extract list with overlapping variants between age at menopause and BCAC data frames
overlap_menopause_bcac = intersect(menopause_data$SNP,
                                   lumA_data$SNP) # 232 SNPs in both data frames


# Extract SNPs that are in age at menarche data frame but NOT in BCAC data frame
proxy_needed_bcac_menopause = setdiff(menopause_data$SNP,
                                      lumA_data$SNP) # 58 SNPs that do not overlap with BCAC

# Check if there are indels among the index variants
table(menopause_data$effect_allele.exposure);
table(menopause_data$other_allele.exposure) # indels snps as index variants; so use selecting_proxies_in_BCAC later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_menopause,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #

# Load list with menopause proxies into R  
proxies_menopause = read.table(file = "combined_query_snp_list.txt",
                               header = T,
                               sep = "\t",
                               row.names = NULL) # n = 64,448

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_menopause = setdiff(proxy_needed_bcac_menopause,
                                 proxies_menopause$query_snp) # 7 SNPs; no duplicates

# Check data
head(proxies_menopause) # correct

table(duplicated(proxies_menopause$query_snp)) # correct; 51 + 7 = 58


# Perform proxy selection for menopause

# Add variables index_SNP and in_BCAC to data
proxies_menopause_2 = variables_proxies(proxies_menopause)

table(proxies_menopause_2$index_snp)
table(proxies_menopause_2$in_BCAC) # 31,460 proxies in BCAC

table(proxies_menopause_2$query_snp==proxies_menopause_2$RS_Number) # correct that there are 44 index snps instead of 51


# Perform actual selection
proxies_menopause_final = selecting_proxies_in_BCAC(proxies_menopause_2) # n = 25 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_menopause_final$query_snp,
                  proxies_menopause_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_menopause_final$RS_Number,
                   menopause_data$SNP) # this is also not the case!

table(duplicated(proxies_menopause_final$RS_Number)) # no duplicates among selected proxies either

rm(check, check2)


# Check R2 values
summary(proxies_menopause_final$R2) # correct, no R2 values below 0.8


# Save menopause proxy dataframe as .txt file
#  write.table(proxies_menopause_final,
#              file = "proxy_search/proxysnps_menopause_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# Clean R environment; remove menopause related objects
rm(menopause_data,
   proxies_menopause,
   proxies_menopause_2,
   proxies_menopause_final,
   overlap_menopause_bcac,
   proxy_needed_bcac_menopause,
   not_searched_menopause)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 11:54h


# 6. Mammographic density (%)

MD_percent_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_percentdensity.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "Beta_PD",
  se_col = "SE_PD",
  eaf_col = "AAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "ref_allele",
  pval_col = "P") # correct, 20 SNPs


# Extract list with overlapping variants between breast density and BCAC data frames
overlap_md_bcac = intersect(MD_percent_data$SNP,
                            lumA_data$SNP) # 12 SNPs in both data frames


# Extract SNPs that are in breast density data frame but NOT in BCAC data frame
proxy_needed_bcac_md = setdiff(MD_percent_data$SNP,
                               lumA_data$SNP) # 8 density SNPs that do not overlap with BCAC


# Check if there are indels among the index variants
table(MD_percent_data$effect_allele.exposure);
table(MD_percent_data$other_allele.exposure) # indels snps as index variants; so use selecting_proxies_in_BCAC later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_md,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #

# Load list with density proxies into R  
proxies_md = read.table(file = "combined_query_snp_list.txt",
                        header = T,
                        sep = "\t",
                        row.names = NULL) # n = 10,726

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_md = setdiff(proxy_needed_bcac_md,
                          proxies_md$query_snp) # 0 SNPs

# Check data
head(proxies_md) # correct

table(duplicated(proxies_md$query_snp)) # correct; 8


# Perform proxy selection for density

# Add variables index_SNP and in_BCAC to data
proxies_md_2 = variables_proxies(proxies_md)

table(proxies_md_2$index_snp)
table(proxies_md_2$in_BCAC) # 8,533 proxies in BCAC


# Perform actual selection
proxies_md_final = selecting_proxies_in_BCAC(proxies_md_2) # n = 7 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_md_final$query_snp,
                  proxies_md_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_md_final$RS_Number,
                   MD_percent_data$SNP) # this is also not the case!


rm(check, check2)


table(duplicated(proxies_md_final$RS_Number)) #there is one duplicated SNP id, which one?

print(proxies_md_final[which(duplicated(proxies_md_final$RS_Number)), c(1:3)]) # rs12048493
print(proxies_md_final[proxies_md_final$RS_Number == "rs12048493", c(1:3)]) # are there other proxies available for query SNP rs11205303?

print(proxies_md[proxies_md$query_snp== "rs11205303" & proxies_md$R2 >= 0.80, c(1:3)]) # only 4 SNPs 
# rs67807996 is the only candidate; R2 is 0.86 and this SNP has not been selected as proxy for another query_snp

# First remove initial row for query_snp rs11205303
proxies_md_final_2 = proxies_md_final %>%
  filter(query_snp != "rs11205303")

# Then, add this the new proxy to proxies_md_final2
proxies_md_final_3 = rbind(proxies_md_final_2,
                           proxies_md[(proxies_md$query_snp == "rs11205303" & proxies_md$RS_Number == "rs67807996"),])

table(duplicated(proxies_md_final_3$RS_Number)) # no duplicates among selected proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_md_final_3$query_snp,
                  proxies_md_final_3$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_md_final_3$RS_Number,
                   MD_percent_data$SNP) # this is also not the case!

rm(check, check2)


# Check R2 values
summary(proxies_md_final_3$R2) # correct, no R2 values below 0.8


# Save density proxy dataframe as .txt file
#  write.table(proxies_md_final_3,
#              file = "proxy_search/proxysnps_density_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # Saved at 20/12/2021


# Clean R environment; remove height related objects
rm(MD_percent_data,
   proxies_md,
   proxies_md_2,
   proxies_md_final,
   proxies_md_final_2,
   proxies_md_final_3,
   overlap_md_bcac,
   proxy_needed_bcac_md,
   not_searched_md)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 13:02h



# 7. Alcohol consumption

alcohol_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_alcohol.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "Alternate Allele Frequency",
  effect_allele_col = "Alternate Allele",
  other_allele_col = "Reference Allele",
  pval_col = "Pvalue",
  samplesize_col = "N") # correct, 99 SNPs


# Extract list with overlapping variants between alcohol consumption and BCAC data frames
overlap_alcohol_bcac = intersect(alcohol_data$SNP,
                                 lumA_data$SNP) # 91 SNPs in both data frames


# Extract SNPs that are in alcohol data frame but NOT in BCAC data frame
proxy_needed_bcac_alcohol = setdiff(alcohol_data$SNP,
                                    lumA_data$SNP) # 8 alcohol SNPs that do not overlap with BCAC

# Drop entry "." to avoid that proxy search gets stuck
proxy_needed_bcac_alcohol_2 = proxy_needed_bcac_alcohol[-8]


# Check if there are indels among the index variants
table(alcohol_data$effect_allele.exposure);
table(alcohol_data$other_allele.exposure) # no indels as index variants; so use selecting_proxies_in_BCAC_excl_indels later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_alcohol_2,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #


# Load list with alcohol proxies into R  
proxies_alcohol = read.table(file = "combined_query_snp_list.txt",
                             header = T,
                             sep = "\t",
                             row.names = NULL) # n = 10,290

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_alcohol = setdiff(proxy_needed_bcac_alcohol,
                               proxies_alcohol$query_snp) # 2 SNPs; no duplicates

# Check data
head(proxies_alcohol) # correct

table(duplicated(proxies_alcohol$query_snp)) # correct; 6 + 2 = 8


# Perform proxy selection for alcohol consumption

# Add variables index_SNP and in_BCAC to data
proxies_alcohol_2 = variables_proxies(proxies_alcohol)

table(proxies_alcohol_2$index_snp)
table(proxies_alcohol_2$in_BCAC) # 5,637 proxies in BCAC

table(proxies_alcohol_2$query_snp==proxies_alcohol_2$RS_Number) # correct that there are 5 index snps instead of 6


# Perform actual selection
proxies_alcohol_final = selecting_proxies_in_BCAC_excl_indels(proxies_alcohol_2) # n = 5 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_alcohol_final$query_snp,
                  proxies_alcohol_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the height_data frame
check2 = intersect(proxies_alcohol_final$RS_Number,
                   alcohol_data$SNP) # this is also not the case!

rm(check, check2)

table(duplicated(proxies_alcohol_final$RS_Number)) #there are no duplicates among the selected proxies


# Check R2 values
summary(proxies_alcohol_final$R2) # correct, no R2 values below 0.8


# Save alcohol proxy dataframe as .txt file
#  write.table(proxies_alcohol_final,
#              file = "proxy_search/proxysnps_alcohol_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# Clean R environment; remove alcohol related objects
rm(alcohol_data,
   proxies_alcohol,
   proxies_alcohol_2,
   proxies_alcohol_final,
   overlap_alcohol_bcac,
   proxy_needed_bcac_alcohol,
   proxy_needed_bcac_alcohol_2,
   not_searched_alcohol)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 13.10h


# 8. Smoking initiation

smoking_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_smoking.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "Alternate Allele Frequency",
  effect_allele_col = "Alternate Allele",
  other_allele_col = "Reference Allele",
  pval_col = "Pvalue") # correct, 378 SNPs


# Extract list with overlapping variants between smoking initiation and BCAC data frames
overlap_smoking_bcac = intersect(smoking_data$SNP,
                                 lumA_data$SNP) # 356 SNPs in both data frames


# Extract SNPs that are in smoking data frame but NOT in BCAC data frame
proxy_needed_bcac_smoking = setdiff(smoking_data$SNP,
                                    lumA_data$SNP) # 22 alcohol SNPs that do not overlap with BCAC

# Check if there are indels among the index variants
table(smoking_data$effect_allele.exposure);
table(smoking_data$other_allele.exposure) # no indels as index variants; so use selecting_proxies_in_BCAC_excl_indels later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_smoking,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #


# Load list with smoking proxies into R  
proxies_smoking = read.table(file = "combined_query_snp_list.txt",
                             header = T,
                             sep = "\t",
                             row.names = NULL) # n = 35,037

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_smoking = setdiff(proxy_needed_bcac_smoking,
                               proxies_smoking$query_snp) # 0 SNPs

# Check data
head(proxies_smoking) # correct

table(duplicated(proxies_smoking$query_snp)) # correct; 22


# Perform proxy selection for smoking

# Add variables index_SNP and in_BCAC to data
proxies_smoking_2 = variables_proxies(proxies_smoking)

table(proxies_smoking_2$index_snp)
table(proxies_smoking_2$in_BCAC) # 19,025 proxies in BCAC

table(proxies_smoking_2$query_snp==proxies_smoking_2$RS_Number) # correct that there are 17 index snps instead of 21


# Perform actual selection
proxies_smoking_final = selecting_proxies_in_BCAC_excl_indels(proxies_smoking_2) # n = 15 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_smoking_final$query_snp,
                  proxies_smoking_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the smoking_data frame
check2 = intersect(proxies_smoking_final$RS_Number,
                   smoking_data$SNP) # this is also not the case!

rm(check, check2)

table(duplicated(proxies_smoking_final$RS_Number)) #there are no duplicates among the selected proxies


# Check R2 values
summary(proxies_smoking_final$R2) # correct, no R2 values below 0.8


# Save smoking proxy dataframe as .txt file
#  write.table(proxies_smoking_final,
#              file = "proxy_search/proxysnps_smoking_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# Clean R environment; remove smoking related objects
rm(smoking_data,
   proxies_smoking,
   proxies_smoking_2,
   proxies_smoking_final,
   overlap_smoking_bcac,
   proxy_needed_bcac_smoking,
   not_searched_smoking)

# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 13:24h


# 9. Physical activity

activity_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/index_stats_activity.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "AF",
  effect_allele_col = "effect_allele",
  other_allele_col = "ref_allele",
  pval_col = "pvalue") # correct, 5 SNPs

# Extract list with overlapping variants between physical activity and BCAC data frames
overlap_activity_bcac = intersect(activity_data$SNP,
                                  lumA_data$SNP) # 3 SNPs in both data frames


# Extract SNPs that are in activity data frame but NOT in BCAC data frame
proxy_needed_bcac_activity = setdiff(activity_data$SNP,
                                     lumA_data$SNP) # 2 activity SNPs that do not overlap with BCAC

# Check if there are indels among the index variants
table(activity_data$effect_allele.exposure);
table(activity_data$other_allele.exposure) # no indels as index variants; so use selecting_proxies_in_BCAC_excl_indels later on!


# Search proxies using LDLink API
LDproxy_batch(snp = proxy_needed_bcac_activity,
              pop = "CEU",
              r2d = "r2",
              token = "1cdbc8b08152",
              append = TRUE) #

# Load list with activity proxies into R  
proxies_activity = read.table(file = "combined_query_snp_list.txt",
                              header = T,
                              sep = "\t",
                              row.names = NULL) # n = 4,101

# How many (and which) SNPs were not available on the 1000G ref panel or a multiallelic variant?
not_searched_activity = setdiff(proxy_needed_bcac_activity,
                                proxies_activity$query_snp) # 0 SNPs

# Check data
head(proxies_activity) # correct

table(duplicated(proxies_activity$query_snp)) # correct; 2


# Perform proxy selection for activity

# Add variables index_SNP and in_BCAC to data
proxies_activity_2 = variables_proxies(proxies_activity)

table(proxies_activity_2$index_snp)
table(proxies_activity_2$in_BCAC) # 3,058 proxies in BCAC

table(proxies_activity_2$query_snp==proxies_activity_2$RS_Number) # correct that there are 0 index snps instead of 2


# Perform actual selection
proxies_activity_final = selecting_proxies_in_BCAC_excl_indels(proxies_activity_2) # n = 2 proxies

# Check if there is overlap between the index snps and searched proxies
check = intersect(proxies_activity_final$query_snp,
                  proxies_activity_final$RS_Number) # no overlap between the two

# Check if there is overlap between the selected proxies and index SNPs in the activity_data frame
check2 = intersect(proxies_activity_final$RS_Number,
                   activity_data$SNP) # this is also not the case!

rm(check, check2)

table(duplicated(proxies_activity_final$RS_Number)) #there are no duplicates among the selected proxies

# Check R2 values
summary(proxies_activity_final$R2) # correct, no R2 values below 0.8


# Save alcohol proxy dataframe as .txt file
#  write.table(proxies_activity_final,
#              file = "proxy_search/proxysnps_activity_bcac_def.txt",
#              quote = F,
#              sep = "\t",
#              row.names = F,
#              col.names = T) # saved at 20/12/2021


# IMPORTANT: also manually remove file "combined_query_snp_list.txt" from folder MR_BCAC before continuing!!
# Note: done at 20/12/2021 at 13:43h

rm(list = ls())

## The end of the script ##
