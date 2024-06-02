## The aim of this Rscript is to make final risk factor data frames for MR ##

library(dplyr)

# Set working directory
setwd("/DATA/users/r.verdiesen/MR_BCAC/")


## load BCAC GWAS data into R ##

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




## Define functions used later on ##

# 1. Function for creating indicator variable proxy_needed
proxy_indicator = function(riskfactor_data) { 
  riskfactor_data %>%
    mutate(proxy_needed = case_when(
      SNP %in% lumA_data$SNP ~ "No",
      ! SNP %in% lumA_data$SNP ~ "Yes"))}


# 2. Function to separate correlated alleles for index and proxy SNPs; using a two-step approach
separate_alleles = function(IVs){
  IVs %>%
    tidyr::separate(Correlated_Alleles, 
                    c("correlated_alleles_A", "correlated_alleles_B"),
                    sep = ",",
                    convert = T,
                    remove = F) %>%
    tidyr::separate(correlated_alleles_A,
                    c("index_allele_A", "proxy_allele_A"),
                    sep = "=",
                    convert = T,
                    remove = F) %>%
    tidyr::separate(correlated_alleles_B,
                    c("index_allele_B", "proxy_allele_B"),
                    sep = "=",
                    convert = T,
                    remove = F)
}


# 3. Function to make indicator variable for palindromic SNPs (index SNPs)
palindromic_indicator_1 = function(IVs){
  IVs %>%
    mutate(index_palindr = case_when(
      ((effect_allele.exposure == "A" & other_allele.exposure == "T") |
         (effect_allele.exposure == "T" & other_allele.exposure == "A") |
         (effect_allele.exposure == "C" & other_allele.exposure == "G") |
         (effect_allele.exposure == "G" & other_allele.exposure == "C")) & 
        (eaf.exposure >= 0.4 & eaf.exposure <= 0.6) ~ 1
    )) 
}


# 4. Function to make indicator variable for palindromic SNPs (proxy SNPs)
palindromic_indicator_2 = function(IVs){
  IVs %>%
    mutate(proxy_palindr = case_when(
      ( Alleles =="(A/T)" | Alleles == "(T/A)" | Alleles == "(C/G)" | Alleles == "(G/C") & (MAF >= 0.4 & MAF <= 0.6) ~ 1
    ))
}


# 5. Function to match alleles
match_alleles = function(IVs){
  IVs %>% 
    mutate(effect_allele.exposure.def = case_when(
      proxy_needed == "Yes" & effect_allele.exposure == index_allele_A ~ proxy_allele_A,
      proxy_needed == "Yes" & effect_allele.exposure == index_allele_B ~ proxy_allele_B,
      proxy_needed == "No" ~ effect_allele.exposure
    )) %>%
    mutate(other_allele.exposure.def = case_when(
      proxy_needed == "Yes" & other_allele.exposure == index_allele_B ~ proxy_allele_B,
      proxy_needed == "Yes" & other_allele.exposure == index_allele_A ~ proxy_allele_A,
      proxy_needed == "No" ~ other_allele.exposure
    ))
} # using this approach, flipping betas is not needed!


# 6. Function to check if matching of alleles went correctly
check_matching = function(matched) {
  matched %>%
    select(SNP, other_allele.exposure, effect_allele.exposure, eaf.exposure,
           proxy_needed, RS_Number, Correlated_Alleles, index_allele_A, proxy_allele_A, index_allele_B, proxy_allele_B,
           effect_allele.exposure.def, other_allele.exposure.def) %>%
    filter(proxy_needed == "Yes")
} # correct


# 7. Function to make final snp id variable
final_rsid = function(matched) {
  matched %>%
    mutate(SNP.def = case_when(
      proxy_needed == "Yes" ~ RS_Number,
      proxy_needed == "No" ~ SNP))
}


# 8. Function to include reason for absence proxy SNP
reason_excl = function(matched){
  matched %>%
    mutate(reason.excl = case_when(
      proxy_needed == "Yes" & is.na(RS_Number) ~ "No proxy SNP available in BCAC subtype summary statistics",
      proxy_needed == "Yes" & index_palindr == 1 ~ "Index SNP palindromic; alleles could not be matched",
      proxy_needed == "Yes" & proxy_palindr == 1 ~ "Proxy SNP palindromic; alleles could not be matched",
      nchar(effect_allele.exposure.def) > 1 | effect_allele.exposure.def == "" | 
        (is.na(effect_allele.exposure.def) & ! is.na(other_allele.exposure.def)) |
        nchar(other_allele.exposure.def) > 1 | other_allele.exposure.def == "" | 
        (is.na(other_allele.exposure.def) & ! is.na(effect_allele.exposure.def)) ~ "Proxy SNP is DEL/INDEL/DELINS, while index SNP is not" # this level should not be present anymore
    ))
}


# 9. Function to create supplemental table 1
supplemental_table1 = function(matched) {
  matched %>%
    select(risk_factor, SNP, proxy_needed, SNP.def, R2, 
           effect_allele.exposure, other_allele.exposure, eaf.exposure,
           effect_allele.exposure.def, other_allele.exposure.def, reason.excl)
}


# 10. Function to create final MR risk factor object
MR_dataframe = function(flipped) {
  flipped %>%
    filter(is.na(reason.excl)) %>% 
    select(SNP.def, effect_allele.exposure.def, other_allele.exposure.def,
           eaf.exposure, samplesize.exposure, beta.exposure,
           se.exposure, pval.exposure, exposure, mr_keep.exposure,
           pval_origin.exposure, id.exposure, data_source.exposure) %>%
    rename(SNP = SNP.def,
           effect_allele.exposure = effect_allele.exposure.def,
           other_allele.exposure = other_allele.exposure.def)
}

# Same function, but for MR dataframes that are lacking a sample size column
MR_dataframe_noN = function(flipped) {
  flipped %>%
    filter(is.na(reason.excl)) %>% 
    select(SNP.def, effect_allele.exposure.def, other_allele.exposure.def,
           eaf.exposure, beta.exposure,
           se.exposure, pval.exposure, exposure, mr_keep.exposure,
           pval_origin.exposure, id.exposure, data_source.exposure) %>%
    rename(SNP = SNP.def,
           effect_allele.exposure = effect_allele.exposure.def,
           other_allele.exposure = other_allele.exposure.def)
}


## Risk factors ##

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


# Load in proxy data
height_proxies = read.table(file = "proxy_search/proxysnps_height_bcac_def.txt",
                            header = T) # correct 270

# Add indicator variable proxy_needed
MR_data = proxy_indicator(height_data)

# Check
table(MR_data$proxy_needed) # correct; 406 proxies needed


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = height_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)
head(MR_incl_proxies) # correct


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies) # no NAs introduced

head(MR_incl_proxies2[MR_incl_proxies2$proxy_needed == "Yes", ]) # correct


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 94 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 12 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # correct
rm(check)


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5)

tail(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Height")

head(MR_incl_proxies7) # correct


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl, 
      useNA = "always") # 3,146 SNPs should remain for analyses

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 35 palindromic index SNPs for which proxies are available

View(palindromic_index[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "proxy_needed",
                           "RS_Number", "Alleles", "MAF", "Correlated_Alleles", "effect_allele.exposure.def", "other_allele.exposure.def",
                           "reason.excl")]) # checked these 35 SNPs manually (via dbSNP); matching went OK for all of them! See Excel file

rm(palindromic_index)


# Create supplemental table 1
height_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(height_supplementtable1) # checked at dd 21.12.2021; correct



# Finalize MR object
height_mrobject_def = MR_dataframe(MR_incl_proxies8) # n = 3,146; correct

head(height_mrobject_def)
table(duplicated(height_mrobject_def$SNP)) # no duplicated SNP ids!


# Remove redundant objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            sure = T)


# 2. BMI

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


# Load in proxy data
bmi_proxies = read.table(file = "proxy_search/proxysnps_bmi_bcac_def.txt",
                         header = T) # correct 85

# Add indicator variable proxy_needed
MR_data = proxy_indicator(bmi_data)

# Check
table(MR_data$proxy_needed) # correct, 106 proxies needed


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = bmi_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)
head(MR_incl_proxies) # correct


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)

head(MR_incl_proxies2[MR_incl_proxies$proxy_needed == "Yes", ]) # correct


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 34 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # 2 palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 0 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # correct
rm(check)


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Body mass index")


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct; 918 SNPs should remain for analysis

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 7 palindromic index SNPs for which proxies are available

View(palindromic_index[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "proxy_needed",
                           "RS_Number", "Alleles", "MAF", "Correlated_Alleles", "effect_allele.exposure.def", "other_allele.exposure.def",
                           "reason.excl")]) # checked these 7 SNPs manually; matching is OK for all of them! See Excel file

rm(palindromic_index)


# Create supplemental table 1
bmi_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(bmi_supplementtable1) # checked at dd 21.12.2021; correct



# Add here step to finalize MR object
bmi_mrobject_def = MR_dataframe(MR_incl_proxies8)# n = 918 = correct

table(duplicated(bmi_mrobject_def$SNP)) # no duplicated SNP ids!


# Remove redundant objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            sure = T)



# 3. T2D

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


# Load in proxy data
t2d_proxies = read.table(file = "proxy_search/proxysnps_t2d_bcac_def.txt",
                         header = T) # correct 29

# Add indicator variable proxy_needed
MR_data = proxy_indicator(t2d_data)

# Check
table(MR_data$proxy_needed) # correct, n = 50


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = t2d_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)

head(MR_incl_proxies) # correct


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)

head(MR_incl_proxies2[MR_incl_proxies2$proxy_needed == "Yes", ]) # correct


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 13 palindromic index SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 0 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # correct
rm(check)


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Type 2 diabetes")

head(MR_incl_proxies7)


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct; 404 SNPs should remain

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 0 palindromic index SNPs for which proxies are available

rm(palindromic_index)


# Create supplemental table 1
t2d_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(t2d_supplementtable1) # checked at dd 21.12.2021; correct



# Finalize MR object
t2d_mrobject_def = MR_dataframe_noN(MR_incl_proxies8) # n = 404; correct

table(duplicated(t2d_mrobject_def$SNP)) # no duplicated SNP ids


# Remove redundant objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def, # still to be made!
            sure = T)


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


# Load in proxy data
menarche_proxies = read.table(file = "proxy_search/proxysnps_menarche_bcac_def.txt",
                              header = T) # correct 24

# Add indicator variable proxy_needed
MR_data = proxy_indicator(menarche_data)

# Check
table(MR_data$proxy_needed) # correct, 44 proxies needed


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = menarche_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)

head(MR_incl_proxies) # correct


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)

head(MR_incl_proxies2[MR_incl_proxies2$proxy_needed == "Yes", ])


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 16 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 3 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # for three INDELS there are missing alleles for the proxy due to the coding; rs376592110, rs10595949, rs375627304
rm(check)

# Subset to index indels for which a proxy is needed
indels_menarche = subset(MR_incl_proxies5,
                         MR_incl_proxies5$proxy_needed == "Yes" &
                           (MR_incl_proxies5$effect_allele.exposure == "D" | MR_incl_proxies5$effect_allele.exposure == "I"))

View(indels_menarche[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", 
                         "RS_Number", "Alleles", "MAF", "Correlated_Alleles")]) # for three of the index indels there are proxies


# Manually change alleles of 3 index INDELS and 1 proxy INDEL
MR_incl_proxies5.1 = MR_incl_proxies5 %>%
  rename(effect_allele.exposure.def.old = effect_allele.exposure.def,
         other_allele.exposure.def.old = other_allele.exposure.def) %>%
  mutate(effect_allele.exposure.def = case_when(
    RS_Number == "rs375627304" | RS_Number == "rs10595949" | RS_Number == "rs376592110" ~ "D",
    RS_Number == "rs200835103" ~ "CG",
    ! RS_Number %in% c("rs375627304", "rs10595949", "rs376592110", "rs200835103") ~ effect_allele.exposure.def.old)) %>%
  mutate(other_allele.exposure.def = case_when(
    RS_Number == "rs375627304" | RS_Number == "rs10595949" | RS_Number == "rs376592110" ~ "I",
    RS_Number == "rs200835103" ~ "G",
    ! RS_Number %in% c("rs375627304", "rs10595949", "rs376592110", "rs200835103") ~ other_allele.exposure.def.old))

# Decision dd 20/10/2021 at 11.33  
# overallbc_data %>% filter(SNP == "rs200835103") # alleles rs200835103 should be changed!

# Check
table(after = MR_incl_proxies5.1$effect_allele.exposure.def,
      before = MR_incl_proxies5.1$effect_allele.exposure.def.old, 
      useNA = "always") # correct

table(after = MR_incl_proxies5.1$other_allele.exposure.def,
      before = MR_incl_proxies5.1$other_allele.exposure.def.old, 
      useNA = "always")

View(MR_incl_proxies5.1[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "proxy_needed", "RS_Number", "Correlated_Alleles", 
                            "effect_allele.exposure.def.old", "other_allele.exposure.def.old",
                            "effect_allele.exposure.def", "other_allele.exposure.def")]) # correct


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5.1)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Age at menarche")

head(MR_incl_proxies7)


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct; 367 SNPs should remain for the analysis

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 6 palindromic index SNPs for which proxies are available

View(palindromic_index[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "proxy_needed",
                           "RS_Number", "Alleles", "MAF", "Correlated_Alleles", "effect_allele.exposure.def", "other_allele.exposure.def",
                           "reason.excl")]) # checked these 6 SNPs manually; matching is OK for all of them! See Excel file

rm(palindromic_index)


# Create supplemental table 1
menarche_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(menarche_supplementtable1) # checked at dd 21.09.2021; correct


# Finalize MR object
menarche_mrobject_def = MR_incl_proxies8 %>%
  filter(reason.excl == "Proxy SNP is DEL/INDEL/DELINS, while index SNP is not" | is.na(reason.excl)) %>% 
  select(SNP.def, effect_allele.exposure.def, other_allele.exposure.def,
         eaf.exposure, samplesize.exposure, beta.exposure,
         se.exposure, pval.exposure, exposure, mr_keep.exposure,
         pval_origin.exposure, id.exposure, data_source.exposure) %>%
  rename(SNP = SNP.def,
         effect_allele.exposure = effect_allele.exposure.def,
         other_allele.exposure = other_allele.exposure.def) # n = 367 = correct; don't forget to manually change supplemental table for rs200835103!!


# Remove redundant objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def,
            menarche_mrobject_def,
            sure = T)


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


# Load in proxy data
menopause_proxies = read.table(file = "proxy_search/proxysnps_menopause_bcac_def.txt",
                               header = T) # correct 25

# Add indicator variable proxy_needed
MR_data = proxy_indicator(menopause_data)

head(MR_data) # correct

# Check
table(MR_data$proxy_needed) # correct; 58


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = menopause_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)

head(MR_incl_proxies) # correct


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)

head(MR_incl_proxies2[MR_incl_proxies2$proxy_needed == "Yes", ]) # correct


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 11 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 0 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # 6 indels (5 index; 1 proxy)
# Index SNPs: rs979481644 (proxy = rs201395656), rs80047570 (proxy = rs34685240), rs375453823 (proxy = rs34696783), rs5800506 (proxy = rs4764695), rs371125724 (proxy = rs71765962)
# Proxy SNPs: rs75066758 (index = rs11245450)
rm(check)

# Check each of these indels in lumA_data to guide manually adjustment below
# bcac = subset(lumA_data,
#              lumA_data$SNP=="rs201395656") # zie aantekening voor handmatig aanpassen op maandag!

# manually adjust alleles for these 6 SNPs
MR_incl_proxies5.1 = MR_incl_proxies5 %>%
  rename(effect_allele.exposure.def.old = effect_allele.exposure.def,
         other_allele.exposure.def.old = other_allele.exposure.def) %>%
  mutate(effect_allele.exposure.def = case_when(
    RS_Number == "rs34685240" ~ "C",
    RS_Number == "rs71765962" ~ "CT",
    RS_Number == "rs75066758" ~ "GAAGAA",
    RS_Number == "rs4764695" ~ "A",
    RS_Number == "rs34696783" ~ "A",
    RS_Number == "rs201395656" ~ "A",
    ! RS_Number %in% c("rs34685240", "rs71765962", "rs75066758", "rs4764695", "rs34696783", "rs201395656") ~ effect_allele.exposure.def.old)) %>%
  mutate(other_allele.exposure.def = case_when(
    RS_Number == "rs34685240" ~ "CT",
    RS_Number == "rs71765962" ~ "C",
    RS_Number == "rs75066758" ~ "G",
    RS_Number == "rs4764695" ~ "G",
    RS_Number == "rs34696783" ~ "AT",
    RS_Number == "rs201395656" ~ "AT",
    ! RS_Number %in% c("rs34685240", "rs71765962", "rs75066758", "rs4764695", "rs34696783", "rs201395656") ~ other_allele.exposure.def.old))

# Check
table(after = MR_incl_proxies5.1$effect_allele.exposure.def,
      before = MR_incl_proxies5.1$effect_allele.exposure.def.old, 
      useNA = "always") # correct

table(after = MR_incl_proxies5.1$other_allele.exposure.def,
      before = MR_incl_proxies5.1$other_allele.exposure.def.old, 
      useNA = "always") # correct

View(MR_incl_proxies5.1[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "proxy_needed", "RS_Number", "Correlated_Alleles", 
                            "effect_allele.exposure.def.old", "other_allele.exposure.def.old",
                            "effect_allele.exposure.def", "other_allele.exposure.def")]) # correct



# Copy effect and other allele variables
MR_incl_proxies5$effect_allele.exposure.def.old = MR_incl_proxies5$effect_allele.exposure.def
MR_incl_proxies5$other_allele.exposure.def.old = MR_incl_proxies5$other_allele.exposure.def


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5.1)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Age at menopause")

head(MR_incl_proxies7)


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct; 257 SNPs can be included

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 3 palindromic index SNPs for which proxies are available

View(palindromic_index[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "proxy_needed",
                           "RS_Number", "Alleles", "MAF", "Correlated_Alleles", "effect_allele.exposure.def", "other_allele.exposure.def",
                           "reason.excl")]) # checked these 3 SNPs manually; matching is OK for all of them! See Excel file

rm(palindromic_index)


# Create supplemental table 1
menopause_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(menopause_supplementtable1) # checked at dd 21.12.2021; correct


# Finalize MR object
menopause_mrobject_def =
  MR_incl_proxies8 %>%
  filter(reason.excl != "No proxy SNP available in BCAC subtype summary statistics" | is.na(reason.excl)) %>% 
  select(SNP.def, effect_allele.exposure.def, other_allele.exposure.def,
         eaf.exposure, beta.exposure,
         se.exposure, pval.exposure, exposure, mr_keep.exposure,
         pval_origin.exposure, id.exposure, data_source.exposure) %>%
  rename(SNP = SNP.def,
         effect_allele.exposure = effect_allele.exposure.def,
         other_allele.exposure = other_allele.exposure.def) # n = 257; correct



# Remove redundant objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def,
            menarche_mrobject_def,
            menopause_mrobject_def,
            sure = T)


# 6. Breast density

density_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_percentdensity.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "Beta_PD",
  se_col = "SE_PD",
  eaf_col = "AAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "ref_allele",
  pval_col = "P") # correct, 20 SNPs


# Load in proxy data
density_proxies = read.table(file = "proxy_search/proxysnps_density_bcac_def.txt",
                             header = T) # correct 7

# Add indicator variable proxy_needed
MR_data = proxy_indicator(density_data)

# Check
table(MR_data$proxy_needed) # correct


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = density_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)

head(MR_incl_proxies2[MR_incl_proxies2$proxy_needed == "Yes", ])


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr);
table(is.na(MR_incl_proxies3$index_palindr)) # 0 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # 1 palindromic proxy SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5) # one INDEL; index SNP = rs34331777 (proxy = rs1561271)

View(check) # correct
rm(check)


# Indel in BCAC to guide manual adjustment below
bcac = subset(lumA_data,
              lumA_data$SNP=="rs1561271") # zie aantekening voor handmatig aanpassen op maandag!

# manually adjust alleles for this SNP
MR_incl_proxies5.1 = MR_incl_proxies5 %>%
  rename(effect_allele.exposure.def.old = effect_allele.exposure.def,
         other_allele.exposure.def.old = other_allele.exposure.def) %>%
  mutate(effect_allele.exposure.def = case_when(
    RS_Number == "rs1561271" ~ "A",
    ! RS_Number %in% c("rs1561271") ~ effect_allele.exposure.def.old)) %>%
  mutate(other_allele.exposure.def = case_when(
    RS_Number == "rs1561271" ~ "G",
    ! RS_Number %in% c("rs1561271") ~ other_allele.exposure.def.old))

# Check
table(after = MR_incl_proxies5.1$effect_allele.exposure.def,
      before = MR_incl_proxies5.1$effect_allele.exposure.def.old, 
      useNA = "always") # correct

table(after = MR_incl_proxies5.1$other_allele.exposure.def,
      before = MR_incl_proxies5.1$other_allele.exposure.def.old, 
      useNA = "always") # correct

View(MR_incl_proxies5.1[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "proxy_needed", "RS_Number", "Correlated_Alleles", 
                            "effect_allele.exposure.def.old", "other_allele.exposure.def.old",
                            "effect_allele.exposure.def", "other_allele.exposure.def")]) # correct


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5.1)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Percent breast density")

head(MR_incl_proxies7)


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct; 18 SNPs to be included in the analyses

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 1 palindromic index SNP for which a proxy is available

View(palindromic_index[, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "proxy_needed",
                           "RS_Number", "Alleles", "MAF", "Correlated_Alleles", "effect_allele.exposure.def", "other_allele.exposure.def",
                           "reason.excl")]) # cannot be checked due to missing eaf 

rm(palindromic_index)


# Create supplemental table 1
density_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(density_supplementtable1) # checked at dd 21.12.2021; correct


# Finalize MR object
density_mrobject_def = MR_incl_proxies8 %>%
  filter(reason.excl == "Proxy SNP is DEL/INDEL/DELINS, while index SNP is not" | is.na(reason.excl)) %>% 
  select(SNP.def, effect_allele.exposure.def, other_allele.exposure.def,
         eaf.exposure, beta.exposure,
         se.exposure, pval.exposure, exposure, mr_keep.exposure,
         pval_origin.exposure, id.exposure, data_source.exposure) %>%
  rename(SNP = SNP.def,
         effect_allele.exposure = effect_allele.exposure.def,
         other_allele.exposure = other_allele.exposure.def) # n = 18; correct

table(duplicated(density_mrobject_def$SNP)) # no duplicated SNPs


# Remove all density related objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def,
            menarche_mrobject_def,
            menopause_mrobject_def,
            density_mrobject_def,
            sure = T)


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


# Load in proxy data
alcohol_proxies = read.table(file = "proxy_search/proxysnps_alcohol_bcac_def.txt",
                             header = T) # correct 5

# Add indicator variable proxy_needed
MR_data = proxy_indicator(alcohol_data)

# Check
table(MR_data$proxy_needed) # correct


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = alcohol_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)

head(MR_incl_proxies) # correct


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)

head(MR_incl_proxies2[MR_incl_proxies2$proxy_needed == "Yes", ])


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 3 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 0 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # correct
rm(check)


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Alcohol consumption")


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 0 palindromic index SNPs for which proxies are available

rm(palindromic_index)


# Create supplemental table 1
alcohol_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(alcohol_supplementtable1) # checked at dd 21.12.2021; correct



# Finalize MR object
alcohol_mrobject_def = MR_dataframe(MR_incl_proxies8) # n = 96; correct


# Remove all alcohol related objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def,
            menarche_mrobject_def,
            menopause_mrobject_def,
            density_mrobject_def,
            alcohol_mrobject_def,
            sure = T)


# 8. Smoking behavior

smoking_data = TwoSampleMR::read_exposure_data(
  filename = "sumstat_riskfactors/cojo_smoking.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "Alternate Allele Frequency",
  effect_allele_col = "Alternate Allele",
  other_allele_col = "Reference Allele",
  pval_col = "Pvalue",
  samplesize_col = "N") # correct, 378 SNPs


# Load in proxy data
smoking_proxies = read.table(file = "proxy_search/proxysnps_smoking_bcac_def.txt",
                             header = T) # correct 15

# Add indicator variable proxy_needed
MR_data = proxy_indicator(smoking_data)

# Check
table(MR_data$proxy_needed) # correct


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = smoking_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 18 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs

# How many of the index SNPs for which proxies are needed are palindromic?
table(MR_incl_proxies4$proxy_needed,
      MR_incl_proxies4$index_palindr) # 0 SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # correct
rm(check)


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5)

head(MR_incl_proxies6[MR_incl_proxies6$proxy_needed == "Yes", c("SNP", "proxy_needed", "RS_Number", "SNP.def")]) # correct

# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Smoking behavior")


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct; 371 SNPs should remain for analysis

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 0 palindromic index SNPs for which proxies are available

rm(palindromic_index)


# Create supplemental table 1
smoking_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(smoking_supplementtable1) # checked at dd 21.12.2021; correct



# Finalize MR object
smoking_mrobject_def = MR_dataframe(MR_incl_proxies8) # n = 371; correct




# Remove all smoking related objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def,
            menarche_mrobject_def,
            menopause_mrobject_def,
            density_mrobject_def,
            alcohol_mrobject_def,
            smoking_mrobject_def,
            sure = T)


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


# Load in proxy data
activity_proxies = read.table(file = "proxy_search/proxysnps_activity_bcac_def.txt",
                              header = T) # correct 273

# Add indicator variable proxy_needed
MR_data = proxy_indicator(activity_data)

# Check
table(MR_data$proxy_needed) # correct


# Merge proxy data to MR_data
MR_incl_proxies = merge(x = MR_data,
                        y = activity_proxies,
                        by.x = "SNP",
                        by.y = "query_snp",
                        all.x = T)


# Separate correlated alleles column into 4 separate columns
MR_incl_proxies2 = separate_alleles(MR_incl_proxies)


# Check how many palindromic SNPs
# Index SNPs
MR_incl_proxies3 = palindromic_indicator_1(MR_incl_proxies2)

table(MR_incl_proxies3$index_palindr) # 0 palindromic SNPs

# Proxy SNPs
MR_incl_proxies4 = palindromic_indicator_2(MR_incl_proxies3)

table(MR_incl_proxies4$proxy_palindr) # suggests that there are no palindromic proxy SNPs, check
table(is.na(MR_incl_proxies4$proxy_palindr)) # indeed no palindromic proxy SNPs


# Match alleles
MR_incl_proxies5 = match_alleles(MR_incl_proxies4)

# Check
check = check_matching(MR_incl_proxies5)

View(check) # correct
rm(check)


# Add column containing final rsid
MR_incl_proxies6 = final_rsid(MR_incl_proxies5)


# Add column containing risk factor
MR_incl_proxies7 = MR_incl_proxies6 %>%
  mutate(risk_factor = "Physical activity")


# Add variable including reason for exclusion before subsequent harmonization of risk factor - outcome data
MR_incl_proxies8 = reason_excl(MR_incl_proxies7)

table(MR_incl_proxies8$reason.excl,
      useNA = "always") # correct, no exclusions have to be made for this phenotype

palindromic_index = subset(MR_incl_proxies8,
                           !is.na(MR_incl_proxies8$RS_Number) &
                             ((MR_incl_proxies8$effect_allele.exposure == "A" & MR_incl_proxies8$other_allele.exposure == "T") | 
                                (MR_incl_proxies8$effect_allele.exposure == "T" & MR_incl_proxies8$other_allele.exposure == "A") |
                                (MR_incl_proxies8$effect_allele.exposure == "C" & MR_incl_proxies8$other_allele.exposure == "G") |
                                (MR_incl_proxies8$effect_allele.exposure == "G" & MR_incl_proxies8$other_allele.exposure == "C"))) 
# 0 palindromic index SNPs for which proxies are available

rm(palindromic_index)


# Create supplemental table 1
activity_supplementtable1 = supplemental_table1(MR_incl_proxies8)

View(activity_supplementtable1) # checked at dd 21.12.2021; correct



# Finalize MR object
activity_mrobject_def = MR_dataframe_noN(MR_incl_proxies8)


# Remove all activity related objects from R environment
gdata::keep(lumA_data,
            check_matching,
            final_rsid,
            match_alleles,
            palindromic_indicator_1,
            palindromic_indicator_2,
            proxy_indicator,
            reason_excl,
            separate_alleles,
            supplemental_table1,
            MR_dataframe,
            MR_dataframe_noN,
            height_mrobject_def,
            bmi_mrobject_def,
            t2d_mrobject_def,
            menarche_mrobject_def,
            menopause_mrobject_def,
            density_mrobject_def,
            alcohol_mrobject_def,
            smoking_mrobject_def,
            activity_mrobject_def,
            sure = T)


# Remove all functions from R environment
rm(list = lsf.str())


# Save curren R environment as .RData object
save.image(file = "MRobjects_subtypespecific_MR_def_version.Rdata") # saved at 21/12/2021


## the end of the script ##
