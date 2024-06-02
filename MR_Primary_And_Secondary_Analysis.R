# This Rscript includes the final two-sample Mendelian randomization analysis for the manuscript "Causal effects of breast cancer risk factors across hormone receptor breast cancer subtypes".
  
#### Exposures of interest:
# Height
# Body mass index
# Type 2 diabetes
# Age at menarche
# Age at menopause
# Percent breast density
# Alcohol consumption
# Smoking behaviour
# Physical activity


#### Outcomes of interest:
# Total breast cancer (all BCAC cases)
# Luminal A-like
# Luminal B/HER2-negative-like
# Luminal B-like
# HER2-enriched-like
# Triple negative


#### Main MR analyses:
# IVW analysis including correlated IVs and LD matrix (figure 1)
  
  
  #### Robust MR analyses (all uncorrelated SNPs):
  # IVW analysis including uncorrelated SNPs only (table 2)
  # MR-Egger analysis (table 2)
  # Weighted median analysis (table 2)
  # Weighted mode analysis (table 2)
  # MR-PRESSO (table 2)
  
  
  #### Sensitivity MR analyses:
  # Analyses including female-specific summary statistics for risk factor (if available) (supplemental table 4)
  # Multivariable MR analyses for BMI & menarche (supplemental figure 3)
  
# Load required libraries
library(dplyr)
library(gdata)
library(ggpubr)
library(ggplot2)
library(inauguration)
library(metafor)
library(LDlinkR)

# First, we load in the required data

# Cleaned summary statistics for the exposures and luminal A-like breast cancer
load("MRobjects_subtypespecific_MR_def_version.Rdata")

# Summary statistics for overall breast cancer and other subtypes
# Overall breast cancer
overallbc_data = TwoSampleMR::read_outcome_data(
   filename = "sumstat_BC/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt",
   sep = " ",
   snp_col = "SNP.Onco",
   beta_col = "Beta.meta",
   se_col = "sdE.meta",
   eaf_col = "EAFcontrols.Onco", # use OncoArray info - because highest nr of participants
   effect_allele_col = "Effect.Meta",
   other_allele_col = "Baseline.Meta",
   pval_col = "p.meta")

# Luminal B/HER2-negative-like
lumBHER2neg_data = TwoSampleMR::read_outcome_data(
   filename = "sumstat_BC/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",
   sep = " ",
   snp_col = "SNP.Onco",
   beta_col = "Luminal_B_HER2Neg_log_or_meta",
   se_col = "Luminal_B_HER2Neg_se_meta",
   eaf_col = "EAFcontrols.Onco", # use OncoArray info - because highest nr of participants
   effect_allele_col = "Effect.Meta",
   other_allele_col = "Baseline.Meta",
   pval_col = "") # no p-value available; TwoSampleMR package infers p-values

# Luminal B-like
lumB_data = TwoSampleMR::read_outcome_data(
   filename = "sumstat_BC/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",
   sep = " ",
   snp_col = "SNP.Onco",
   beta_col = "Luminal_B_log_or_meta",
   se_col = "Luminal_B_se_meta",
   eaf_col = "EAFcontrols.Onco", # use OncoArray info - because highest nr of participants
   effect_allele_col = "Effect.Meta",
   other_allele_col = "Baseline.Meta",
   pval_col = "") # no p-value available; TwoSampleMR package infers p-values

# HER2-enriched-like
HER2enr_data = TwoSampleMR::read_outcome_data(
   filename = "sumstat_BC/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",
   sep = " ",
   snp_col = "SNP.Onco",
   beta_col = "HER2_Enriched_log_or_meta",
   se_col = "HER2_Enriched_se_meta",
   eaf_col = "EAFcontrols.Onco", # use OncoArray info - because highest nr of participants
   effect_allele_col = "Effect.Meta",
   other_allele_col = "Baseline.Meta",
   pval_col = "") # no p-value available; TwoSampleMR package infers p-values

# Triple negative
TNBC_data = TwoSampleMR::read_outcome_data(
   filename = "sumstat_BC/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",
   sep = " ",
   snp_col = "SNP.Onco",
   beta_col = "Triple_Neg_log_or_meta",
   se_col = "Triple_Neg_se_meta",
   eaf_col = "EAFcontrols.Onco", # use OncoArray info - because highest nr of participants
   effect_allele_col = "Effect.Meta",
   other_allele_col = "Baseline.Meta",
   pval_col = "") # no p-value available; TwoSampleMR package infers p-values


# The format of the SNP column for overall breast cancer, luminal B/HER2-negative-like, luminal B-like, HER2-enriched, and TNBC is incorrect; it's in the rsid:bp:A1:A2 format.
 
# Therefore we first have to create a new SNP variable including only the rsid,
# for this we use the sapply:

# Overall breast cancer
overallbc_data$SNP_old = overallbc_data$SNP

overallbc_data$SNP = sapply(strsplit(overallbc_data$SNP_old, 
                                     split = ":"), "[", 1) # this does not result in an rsid for every entry;

head(overallbc_data) # correct
overallbc_data$SNP_old = NULL


# Luminal B/HER2-negative-like
lumBHER2neg_data$SNP_old = lumBHER2neg_data$SNP

lumBHER2neg_data$SNP = sapply(strsplit(lumBHER2neg_data$SNP_old, 
                                       split = ":"), "[", 1) 

head(lumBHER2neg_data) # correct
lumBHER2neg_data$SNP_old = NULL


# Luminal B-like
lumB_data$SNP_old = lumB_data$SNP

lumB_data$SNP = sapply(strsplit(lumB_data$SNP_old, 
                                split = ":"), "[", 1)

head(lumB_data) # correct
lumB_data$SNP_old = NULL


# HER2-enriched
HER2enr_data$SNP_old = HER2enr_data$SNP

HER2enr_data$SNP = sapply(strsplit(HER2enr_data$SNP_old, 
                                   split = ":"), "[", 1)

head(HER2enr_data) # correct
HER2enr_data$SNP_old = NULL


# TNBC
TNBC_data$SNP_old = TNBC_data$SNP

TNBC_data$SNP = sapply(strsplit(TNBC_data$SNP_old, 
                                split = ":"), "[", 1)

head(TNBC_data) # correct
TNBC_data$SNP_old = NULL


# Save objects generated with code above in manuscript folder
save.image(file = "manuscript/allMRdata_finalversion.Rdata") 

# Load MR sumstats into R before proceeding
load(file = "manuscript/allMRdata_finalversion.Rdata")


################################################################################
# Create a first MR pipeline that:
# 1. Harmonizes the exposure-outcome data
# 2. Converts the harmonized data to an MRinput object (includes calculation LD matrix)
# 3. Performs a IVW analysis with multiplicative effects including a LD matrix
# 4. Saves the results from these analyses

mr_pipeline_ld <- function(riskfactor_data, outcome_data, riskfactor_label, outcome_label, analysis) {
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = riskfactor_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  # Exclude palindromic SNPs with intermediate allele frequencies before converting to MRinput
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Conversion to MRinput object including calculation of LD matrix
  mrinput <- TwoSampleMR::dat_to_MRInput(dat = harmon_data_2, #
                                         get_correlations = TRUE, # calculate LD matrix
                                         pop = "EUR")
  
  # Extract object required for MR analyses from list "mrinput"
  mr_data = mrinput[[1]]
  
  
  # Create a matrix to save the results
  results = as.data.frame(matrix(nrow = 1,
                          ncol = 9,
                          dimnames = list(1, 
                                          c("Risk factor", "Outcome", "n.SNPs", "IVW.est", "se.est", "CI.low", "CI.high", "p.hetr.snps", "Analysis"))))
  
  results[1,1] = riskfactor_label
  
  results[1,2] = outcome_label
  
  
  # Perform IVW analysis with multiplicative effects
  ivw_ld <- MendelianRandomization::mr_ivw(mr_data)
  
    # Save results in matrix
    results[1,3] = ivw_ld@SNPs
  
    results[1,4] = ivw_ld@Estimate
    
    results[1,5] = ivw_ld@StdError
                          
    results[1,6] = ivw_ld@CILower

    results[1,7] = ivw_ld@CIUpper
    
    results[1,8] = round(ivw_ld@Heter.Stat[2], digits = 3)
    
    results[1,9] = analysis
    
    return(results)
  
}


mr_pipeline_pca <- function(riskfactor_data, outcome_data, riskfactor_label, outcome_label, analysis) {
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = riskfactor_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  # Exclude palindromic SNPs with intermediate allele frequencies before converting to MRinput
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Calculate LD matrix
  rho = ieugwasr::ld_matrix(variants = harmon_data_2$SNP,
                          with_alleles = FALSE,
                          plink_bin = genetics.binaRies::get_plink_binary(),
                          bfile = "1000G_ref/EUR")
  
  # Subset data to SNPs in correlation matrix
  harmon_data_3 <- harmon_data_2 %>%
  filter(SNP %in% row.names(rho))
  
  
  # Define vectors for PCA IVW MR
  betaXG = harmon_data_3$beta.exposure
  
  sebetaYG = harmon_data_3$se.outcome

  betaYG = harmon_data_3$beta.outcome
  
  # Create a matrix to save the results
  results = as.data.frame(matrix(nrow = 1,
                          ncol = 9,
                          dimnames = list(1, 
                                          c("Risk factor", "Outcome", "n.SNPs", "IVW.est", "se.est", "CI.low", "CI.high", "p.hetr.snps", "Analysis"))))
  
  results[1,1] = riskfactor_label
  
  results[1,2] = outcome_label
  
  
  # Calculate IVW estimate (accounting for correlation) using principal components:

  Phi = (betaXG/sebetaYG)%o%(betaXG/sebetaYG)*rho

    summary(prcomp(Phi, scale=FALSE))

  K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]
  # K is number of principal components to include in analysis

  # The code below includes principal components to explain 99% of variance in the risk factor

  betaXG0 = as.numeric(betaXG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])

  betaYG0 = as.numeric(betaYG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])

  Omega = sebetaYG%o%sebetaYG*rho

  pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]

  beta_IVWcorrel.pc = solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0

  se_IVWcorrel.fixed.pc = sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0))
  
  resid = betaYG0-beta_IVWcorrel.pc*betaXG0
  # checked warning; can be ignored same values are produced when both beta_IVWcorrel.pc and betaXG0 are explicitly specified as vectors
    
  se_IVWcorrel.random.pc = sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0))*max(sqrt(t(resid)%*%solve(pcOmega)%*%resid/(length(betaXG0)-1)),1)
  
  
  # Save results
  results[1,3] = K # number of included PCs instead of number of SNPs!
  
  results[1,4] = beta_IVWcorrel.pc
  
  results[1,5] = se_IVWcorrel.random.pc 
  # This is the difference in file from dd02022022 and dd18012022; dd02022022 version contains random se's; dd18012022 version contains fixed se's
                          
  results[1,6] = beta_IVWcorrel.pc - (se_IVWcorrel.random.pc*1.96) 
  
  results[1,7] = beta_IVWcorrel.pc + (se_IVWcorrel.random.pc*1.96)
  
  results[1,9] = analysis
  
  
  return(results)
  
}
  

################################################################################  
# Create a third MR pipeline that:
# 1. Clumps the genetic instruments
# 2. Harmonizes the clumped exposure-outcome data
# 3. Converts the harmonized data to an MRinput object
# 4. Performs an IVW analysis with multiplicative effects
# 5. Performs a MR-Egger analysis
# 6. Performs a weighted median analysis
# 7. Performs a weighted mode analysis
# 8. Performs a MR-PRESSO analysis
# 9. Saves the results from these analyses

mr_pipeline_clumped <- function(riskfactor_data, outcome_data, riskfactor_label, outcome_label){
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = riskfactor_data,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before converting to MRinput
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Conversion to MRinput object
  mrinput <- TwoSampleMR::dat_to_MRInput(dat = harmon_data_2, #
                                         get_correlations = FALSE) # no need to calculate LD matrix
  
  # Extract object required for MR analyses from list "mrinput"
  mr_data = mrinput[[1]]
  
  
  # Create a matrix to save the results
  results = as.data.frame(matrix(nrow = 1,
                          ncol = 23,
                          dimnames = list(1, 
                                          c("Risk factor", "Outcome", "n.SNPs", "IVW.est", "IVW.se", "CI.low", "CI.high", "p.hetr.snps",
                                            "Egger.est", "Egger.se", "Egger.CI", "P.pleiotr", "Median.est", "Median.se", "Median.CI",
                                            "Mode.est", "Mode.se", "Mode.CI", "PRESSO.est", "PRESSO.se", "PRESSO.CI", 
                                            "PRESSO.pleio.p", "PRESSO.n.outliers"))))
  
  results[1,1] = riskfactor_label
  
  results[1,2] = outcome_label
  
  
  # Perform IVW analysis with multiplicative effects
  ivw_clumped <- MendelianRandomization::mr_ivw(mr_data)
  
    # Save results in matrix
    results[1,3] = ivw_clumped@SNPs # these we convert later to ORs for table 2
  
    results[1,4] = ivw_clumped@Estimate
    
    results[1,5] = ivw_clumped@StdError
                          
    results[1,6] = ivw_clumped@CILower

    results[1,7] = ivw_clumped@CIUpper
    
    results[1,8] = round(ivw_clumped@Heter.Stat[2], digits = 3)
    
    
  # Perform MR-Egger
  egger_clumped <- MendelianRandomization::mr_egger(mr_data)
  
    # Save results in matrix
    results[1,9] = egger_clumped@Estimate
    
    results[1,10] = egger_clumped@StdError.Est
    
    results[1,11] = paste(egger_clumped@CILower.Est, egger_clumped@CIUpper.Est, sep = ", ")

    results[1,12] = round(egger_clumped@Pleio.pval, digits = 3)
    
    
  # Perform weighted median
  median_clumped <- MendelianRandomization::mr_median(mr_data)
  
    # Save results in matrix
    results[1,13] = median_clumped@Estimate
    
    results[1,14] = median_clumped@StdError
                      
    results[1,15] = paste(median_clumped@CILower, median_clumped@CIUpper, sep = ", ")
    
    
  # Perform weighted mode
  mode_clumped <- TwoSampleMR::mr_mode(dat = harmon_data_2, # attention: through TwoSampleMR package; so use data that is harmonized data before conversion!
                                       mode_method = "Weighted mode")

    # Save results in matrix
    results[1,16] = mode_clumped$b
    
    results[1,17] = mode_clumped$se
    
    results[1,18] = paste(mode_clumped$b - (mode_clumped$se*1.96), 
                          mode_clumped$b + (mode_clumped$se*1.96), sep = ", ")
    
    
  # Perform MR-PRESSO via R base in terminal! See separate scripts for that; was needed to run analyses in parallel, otherwise MR-PRESSO would have taken too long! 
   
    return(results)
    
}


mr_pipeline_clumped_female <- function(clumped_data, outcome_data, riskfactor_label, outcome_label){
  
  # Clumping of genetic instruments for risk factor is done outside of this function; same SNPs are included as after clumping of sex-combined data!
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before converting to MRinput
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Conversion to MRinput object
  mrinput <- TwoSampleMR::dat_to_MRInput(dat = harmon_data_2, #
                                         get_correlations = FALSE) # no need to calculate LD matrix
  
  
  # Extract object required for MR analyses from list "mrinput"
  mr_data = mrinput[[1]]
  
  
  # Create a matrix to save the results
  results = as.data.frame(matrix(nrow = 1,
                                 ncol = 23,
                                 dimnames = list(1, 
                                                 c("Risk factor", "Outcome", "n.SNPs", "IVW.est", "IVW.se", "CI.low", "CI.high", "p.hetr.snps",
                                                   "Egger.est", "Egger.se", "Egger.CI", "P.pleiotr", "Median.est", "Median.se", "Median.CI",
                                                   "Mode.est", "Mode.se", "Mode.CI", "PRESSO.est", "PRESSO.se", "PRESSO.CI", 
                                                   "PRESSO.pleio.p", "PRESSO.n.outliers"))))
  
  results[1,1] = riskfactor_label
  
  results[1,2] = outcome_label
  
  
  # Perform IVW analysis with multiplicative effects
  ivw_clumped <- MendelianRandomization::mr_ivw(mr_data)
  
  # Save results in matrix
  results[1,3] = ivw_clumped@SNPs # these we convert later to ORs for table 2
  
  results[1,4] = ivw_clumped@Estimate
  
  results[1,5] = ivw_clumped@StdError
  
  results[1,6] = ivw_clumped@CILower
  
  results[1,7] = ivw_clumped@CIUpper
  
  results[1,8] = round(ivw_clumped@Heter.Stat[2], digits = 3)
  
  
  # Perform MR-Egger
  egger_clumped <- MendelianRandomization::mr_egger(mr_data)
  
  # Save results in matrix
  results[1,9] = egger_clumped@Estimate
  
  results[1,10] = egger_clumped@StdError.Est
  
  results[1,11] = paste(egger_clumped@CILower.Est, egger_clumped@CIUpper.Est, sep = ", ")
  
  results[1,12] = round(egger_clumped@Pleio.pval, digits = 3)
  
  
  # Perform weighted median
  median_clumped <- MendelianRandomization::mr_median(mr_data)
  
  # Save results in matrix
  results[1,13] = median_clumped@Estimate
  
  results[1,14] = median_clumped@StdError
  
  results[1,15] = paste(median_clumped@CILower, median_clumped@CIUpper, sep = ", ")
  
  
  tryCatch({
  
  # Perform weighted mode
  mode_clumped <- TwoSampleMR::mr_mode(dat = harmon_data_2, # attention: through TwoSampleMR package; so use data that is harmonized data before conversion!
                                       mode_method = "Weighted mode")
  
  # Save results in matrix
  results[1,16] = mode_clumped$b
  
  results[1,17] = mode_clumped$se
  
  results[1,18] = paste(mode_clumped$b - (mode_clumped$se*1.96), 
                        mode_clumped$b + (mode_clumped$se*1.96), sep = ", ")
  
  },
  error = function(e) {
    message('An Error Occurred')
    print(e)
    return(results)
  })}


  # Perform MR-PRESSO via R base in terminal!
  
  
  return(results)
  


#################################### height ####################################

# Run main analyses for height; IVW analysis including LD matrix not possible due to correlation SNPs
height_overall_pca <- mr_pipeline_pca(riskfactor_data = height_mrobject_def, 
                                      outcome_data = overallbc_data, 
                                      riskfactor_label = "Height", 
                                      outcome_label = "All BCAC breast cancer cases",
                                      analysis = "IVW including PCs")

height_lumA_pca <- mr_pipeline_pca(riskfactor_data = height_mrobject_def, 
                                   outcome_data = lumA_data, 
                                   riskfactor_label = "Height", 
                                   outcome_label = "Luminal A-like",
                                   analysis = "IVW including PCs")

height_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = height_mrobject_def, 
                                          outcome_data = lumBHER2neg_data, 
                                          riskfactor_label = "Height", 
                                          outcome_label = "Luminal B/HER2-negative like",
                                          analysis = "IVW including PCs")

height_lumB_pca <- mr_pipeline_pca(riskfactor_data = height_mrobject_def, 
                                   outcome_data = lumB_data, 
                                   riskfactor_label = "Height", 
                                   outcome_label = "Luminal B-like",
                                   analysis = "IVW including PCs")

height_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = height_mrobject_def, 
                                      outcome_data = HER2enr_data, 
                                      riskfactor_label = "Height", 
                                      outcome_label = "HER2-enriched-like",
                                      analysis = "IVW including PCs")

height_TNBC_pca <- mr_pipeline_pca(riskfactor_data = height_mrobject_def, 
                                   outcome_data = TNBC_data, 
                                   riskfactor_label = "Height", 
                                   outcome_label = "Triple negative",
                                   analysis = "IVW including PCs")


height_all_fig1 <- rbind(height_overall_pca,
                         height_lumA_pca,
                         height_lumBHER2neg_pca,
                         height_lumB_pca,
                         height_HER2enr_pca,
                         height_TNBC_pca)

height_all_fig1

rm(height_overall_pca, height_lumA_pca, height_lumBHER2neg_pca, height_lumB_pca, height_HER2enr_pca, height_TNBC_pca)

# Save height_all_fig1 in separate Rdata object
# save(height_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Run pipeline for height with clumped variants 
height_overall_clumped <- mr_pipeline_clumped(riskfactor_data = height_mrobject_def, 
                                              outcome_data = overallbc_data, 
                                              riskfactor_label = "Height", 
                                              outcome_label = "All BCAC breast cancer cases")

height_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = height_mrobject_def, 
                                           outcome_data = lumA_data, 
                                           riskfactor_label = "Height", 
                                           outcome_label = "Luminal A-like")

height_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = height_mrobject_def, 
                                                  outcome_data = lumBHER2neg_data, 
                                                  riskfactor_label = "Height", 
                                                  outcome_label = "Luminal B/HER2-negative like")

height_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = height_mrobject_def, 
                                           outcome_data = lumB_data,
                                           riskfactor_label = "Height", 
                                           outcome_label = "Luminal B-like")

height_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = height_mrobject_def, 
                                              outcome_data = HER2enr_data, 
                                              riskfactor_label = "Height", 
                                              outcome_label = "HER2-enriched-like")

height_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = height_mrobject_def, 
                                           outcome_data = TNBC_data, 
                                           riskfactor_label = "Height", 
                                           outcome_label = "Triple negative")

height_all_table2 <- rbind(height_overall_clumped,
                           height_lumA_clumped,
                           height_lumBHER2neg_clumped,
                           height_lumB_clumped,
                           height_HER2enr_clumped,
                           height_TNBC_clumped)

height_all_table2

rm(height_overall_clumped, height_lumA_clumped, height_lumBHER2neg_clumped, height_lumB_clumped, height_HER2enr_clumped, height_TNBC_clumped)



# Run primary analyses for height using female-specific weights for genetic instruments

# Load female specific MR object into R
load("sumstat_riskfactors/female_betas/height_mrobject_women.RData")


height_overall_pca_female <- mr_pipeline_pca(riskfactor_data = height_mrobject_women, 
                                             outcome_data = overallbc_data, 
                                             riskfactor_label = "Height", 
                                             outcome_label = "All BCAC breast cancer cases",
                                             analysis = "IVW including PCs; female betas")

height_lumA_pca_female <- mr_pipeline_pca(riskfactor_data = height_mrobject_women, 
                                          outcome_data = lumA_data, 
                                          riskfactor_label = "Height", 
                                          outcome_label = "Luminal A-like",
                                          analysis = "IVW including PCs; female betas")

height_lumBHER2neg_pca_female <- mr_pipeline_pca(riskfactor_data = height_mrobject_women, 
                                                 outcome_data = lumBHER2neg_data, 
                                                 riskfactor_label = "Height", 
                                                 outcome_label = "Luminal B/HER2-negative like",
                                                 analysis = "IVW including PCs; female betas")

height_lumB_pca_female <- mr_pipeline_pca(riskfactor_data = height_mrobject_women, 
                                          outcome_data = lumB_data, 
                                          riskfactor_label = "Height", 
                                          outcome_label = "Luminal B-like",
                                          analysis = "IVW including PCs; female betas")

height_HER2enr_pca_female <- mr_pipeline_pca(riskfactor_data = height_mrobject_women, 
                                             outcome_data = HER2enr_data, 
                                             riskfactor_label = "Height", 
                                             outcome_label = "HER2-enriched-like",
                                             analysis = "IVW including PCs; female betas")

height_TNBC_pca_female <- mr_pipeline_pca(riskfactor_data = height_mrobject_women, 
                                          outcome_data = TNBC_data, 
                                          riskfactor_label = "Height", 
                                          outcome_label = "Triple negative",
                                          analysis = "IVW including PCs; female betas")


height_pca_all_female <- rbind(height_overall_pca_female,
                               height_lumA_pca_female,
                               height_lumBHER2neg_pca_female,
                               height_lumB_pca_female,
                               height_HER2enr_pca_female,
                               height_TNBC_pca_female)

height_pca_all_female

rm(height_overall_pca_female, height_lumA_pca_female, height_lumBHER2neg_pca_female, height_lumB_pca_female, height_HER2enr_pca_female, height_TNBC_pca_female)


# Save height_pca_all_female in separate Rdata object
# save(height_pca_all_female, file = "manuscript/data_height_femalespecific.Rdata") # saved 01-12-2022 at 11:56

rm(height_pca_all_female)


# Repeat analyses with clumped variants with female-specific estimates
# Maybe best to do clumping based on sex-combined data; use that dataframe to subset female-specific dataframe
# and subsequently run the pipeline without clumping; because that has then already been "done"

height_mrobject_clumped <- TwoSampleMR::clump_data(dat = height_mrobject_def,
                                                   clump_kb = 1000,
                                                   clump_r2 = 0.001,
                                                   clump_p1 = 4.99e-08,
                                                   clump_p2 = 4.99e-08,
                                                   pop = "EUR")


height_mrobject_women_clumped <- height_mrobject_women %>%
  filter(SNP %in% height_mrobject_clumped$SNP)

rm(height_mrobject_clumped,
   height_mrobject_def,
   height_mrobject_women)


height_overall_clumped_female <- mr_pipeline_clumped_female(clumped_data = height_mrobject_women_clumped, 
                                                            outcome_data = overallbc_data, 
                                                            riskfactor_label = "Height", 
                                                            outcome_label = "All BCAC breast cancer cases")

height_lumA_clumped_female <- mr_pipeline_clumped_female(clumped_data = height_mrobject_women_clumped, 
                                                         outcome_data = lumA_data, 
                                                         riskfactor_label = "Height", 
                                                         outcome_label = "Luminal A-like")

height_lumBHER2neg_clumped_female <- mr_pipeline_clumped_female(clumped_data = height_mrobject_women_clumped, 
                                                                outcome_data = lumBHER2neg_data, 
                                                                riskfactor_label = "Height", 
                                                                outcome_label = "Luminal B/HER2-negative like")

height_lumB_clumped_female <- mr_pipeline_clumped_female(clumped_data = height_mrobject_women_clumped, 
                                                         outcome_data = lumB_data,
                                                         riskfactor_label = "Height", 
                                                         outcome_label = "Luminal B-like")

height_HER2enr_clumped_female <- mr_pipeline_clumped_female(clumped_data = height_mrobject_women_clumped, 
                                                            outcome_data = HER2enr_data, 
                                                            riskfactor_label = "Height", 
                                                            outcome_label = "HER2-enriched-like")

height_TNBC_clumped_female <- mr_pipeline_clumped_female(clumped_data = height_mrobject_women_clumped, 
                                                         outcome_data = TNBC_data, 
                                                         riskfactor_label = "Height", 
                                                         outcome_label = "Triple negative")

height_all_clumped_female <- rbind(height_overall_clumped_female,
                                   height_lumA_clumped_female,
                                   height_lumBHER2neg_clumped_female,
                                   height_lumB_clumped_female,
                                   height_HER2enr_clumped_female,
                                   height_TNBC_clumped_female)

height_all_clumped_female

rm(height_overall_clumped_female, height_lumA_clumped_female, height_lumBHER2neg_clumped_female, height_lumB_clumped_female, height_HER2enr_clumped_female, height_TNBC_clumped_female)

rm(height_mrobject_women_clumped)


# Save height_all_clumped_female in separate Rdata object
# save(height_all_clumped_female, file = "manuscript/data_height_femalespecific_clumped.Rdata") # saved dd 01.12.2022 14:38

rm(height_all_clumped_female)



##################################### BMI ######################################

# Run main analyses for BMI; IVW analysis including LD matrix not possible due to correlation SNPs
bmi_overall_pca <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_def, 
                                   outcome_data = overallbc_data, 
                                   riskfactor_label = "Body mass index (kg/m2)", 
                                   outcome_label = "All BCAC breast cancer cases",
                                   analysis = "IVW including PCs")

bmi_lumA_pca <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_def, 
                                outcome_data = lumA_data, 
                                riskfactor_label = "Body mass index (kg/m2)", 
                                outcome_label = "Luminal A-like",
                                analysis = "IVW including PCs")

bmi_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_def, 
                                       outcome_data = lumBHER2neg_data, 
                                       riskfactor_label = "Body mass index (kg/m2)", 
                                       outcome_label = "Luminal B/HER2-negative like",
                                       analysis = "IVW including PCs")

bmi_lumB_pca <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_def, 
                                outcome_data = lumB_data, 
                                riskfactor_label = "Body mass index (kg/m2)", 
                                outcome_label = "Luminal B-like",
                                analysis = "IVW including PCs")

bmi_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_def, 
                                   outcome_data = HER2enr_data, 
                                   riskfactor_label = "Body mass index (kg/m2)", 
                                   outcome_label = "HER2-enriched-like",
                                   analysis = "IVW including PCs")

bmi_TNBC_pca <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_def, 
                                outcome_data = TNBC_data, 
                                riskfactor_label = "Body mass index (kg/m2)", 
                                outcome_label = "Triple negative",
                                analysis = "IVW including PCs")


bmi_all_fig1 <- rbind(bmi_overall_pca,
                      bmi_lumA_pca,
                      bmi_lumBHER2neg_pca,
                      bmi_lumB_pca,
                      bmi_HER2enr_pca,
                      bmi_TNBC_pca)

bmi_all_fig1

rm(bmi_overall_pca, bmi_lumA_pca, bmi_lumBHER2neg_pca, bmi_lumB_pca, bmi_HER2enr_pca, bmi_TNBC_pca)

# Save bmi_all_fig1 in separate Rdata object
# cgwtools::resave(bmi_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Run pipeline for BMI with clumped variants
bmi_overall_clumped <- mr_pipeline_clumped(riskfactor_data = bmi_mrobject_def, 
                                           outcome_data = overallbc_data, 
                                           riskfactor_label = "Body mass index (kg/m2)", 
                                           outcome_label = "All BCAC breast cancer cases")

bmi_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = bmi_mrobject_def, 
                                        outcome_data = lumA_data, 
                                        riskfactor_label = "Body mass index (kg/m2)", 
                                        outcome_label = "Luminal A-like")

bmi_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = bmi_mrobject_def, 
                                               outcome_data = lumBHER2neg_data, 
                                               riskfactor_label = "Body mass index (kg/m2)", 
                                               outcome_label = "Luminal B/HER2-negative like")

bmi_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = bmi_mrobject_def, 
                                        outcome_data = lumB_data,
                                        riskfactor_label = "Body mass index (kg/m2)", 
                                        outcome_label = "Luminal B-like")

bmi_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = bmi_mrobject_def, 
                                           outcome_data = HER2enr_data, 
                                           riskfactor_label = "Body mass index (kg/m2)", 
                                           outcome_label = "HER2-enriched-like")

bmi_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = bmi_mrobject_def, 
                                        outcome_data = TNBC_data, 
                                        riskfactor_label = "Body mass index (kg/m2)", 
                                        outcome_label = "Triple negative")

bmi_all_table2 <- rbind(bmi_overall_clumped,
                        bmi_lumA_clumped,
                        bmi_lumBHER2neg_clumped,
                        bmi_lumB_clumped,
                        bmi_HER2enr_clumped,
                        bmi_TNBC_clumped)

bmi_all_table2

rm(bmi_overall_clumped, bmi_lumA_clumped, bmi_lumBHER2neg_clumped, bmi_lumB_clumped, bmi_HER2enr_clumped, bmi_TNBC_clumped)


# Run primary analyses for bmi using female-specific weights for genetic instruments

# Load female specific MR object into R
load("sumstat_riskfactors/female_betas/bmi_mrobject_women.RData")


bmi_overall_pca_female <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_women, 
                                          outcome_data = overallbc_data, 
                                          riskfactor_label = "Body mass index (kg/m2)", 
                                          outcome_label = "All BCAC breast cancer cases",
                                          analysis = "IVW including PCs; female betas")

bmi_lumA_pca_female <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_women, 
                                       outcome_data = lumA_data, 
                                       riskfactor_label = "Body mass index (kg/m2)", 
                                       outcome_label = "Luminal A-like",
                                       analysis = "IVW including PCs; female betas")

bmi_lumBHER2neg_pca_female <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_women, 
                                              outcome_data = lumBHER2neg_data, 
                                              riskfactor_label = "Body mass index (kg/m2)", 
                                              outcome_label = "Luminal B/HER2-negative like",
                                              analysis = "IVW including PCs; female betas")

bmi_lumB_pca_female <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_women, 
                                       outcome_data = lumB_data, 
                                       riskfactor_label = "Body mass index (kg/m2)", 
                                       outcome_label = "Luminal B-like",
                                       analysis = "IVW including PCs; female betas")

bmi_HER2enr_pca_female <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_women, 
                                          outcome_data = HER2enr_data, 
                                          riskfactor_label = "Body mass index (kg/m2)", 
                                          outcome_label = "HER2-enriched-like",
                                          analysis = "IVW including PCs; female betas")

bmi_TNBC_pca_female <- mr_pipeline_pca(riskfactor_data = bmi_mrobject_women, 
                                       outcome_data = TNBC_data, 
                                       riskfactor_label = "Body mass index (kg/m2)", 
                                       outcome_label = "Triple negative",
                                       analysis = "IVW including PCs; female betas")


bmi_pca_all_female <- rbind(bmi_overall_pca_female,
                            bmi_lumA_pca_female,
                            bmi_lumBHER2neg_pca_female,
                            bmi_lumB_pca_female,
                            bmi_HER2enr_pca_female,
                            bmi_TNBC_pca_female)

bmi_pca_all_female

rm(bmi_overall_pca_female, bmi_lumA_pca_female, bmi_lumBHER2neg_pca_female, bmi_lumB_pca_female, bmi_HER2enr_pca_female, bmi_TNBC_pca_female)


# Save bmi_pca_all_female in separate Rdata object
# save(bmi_pca_all_female, file = "manuscript/data_bmi_femalespecific.Rdata") # saved 01-12-2022 at 15:03

rm(bmi_pca_all_female)


# Repeat analyses with clumped variants with female-specific estimates
# Maybe best to do clumping based on sex-combined data; use that dataframe to subset female-specific dataframe
# and subsequently run the pipeline without clumping; because that has then already been "done"

bmi_mrobject_clumped <- TwoSampleMR::clump_data(dat = bmi_mrobject_def,
                                                clump_kb = 1000,
                                                clump_r2 = 0.001,
                                                clump_p1 = 4.99e-08,
                                                clump_p2 = 4.99e-08,
                                                pop = "EUR")


bmi_mrobject_women_clumped <- bmi_mrobject_women %>%
  filter(SNP %in% bmi_mrobject_clumped$SNP)

rm(bmi_mrobject_clumped,
   bmi_mrobject_def,
   bmi_mrobject_women)


bmi_overall_clumped_female <- mr_pipeline_clumped_female(clumped_data = bmi_mrobject_women_clumped, 
                                                         outcome_data = overallbc_data, 
                                                         riskfactor_label = "Body mass index (kg/m2)", 
                                                         outcome_label = "All BCAC breast cancer cases")

bmi_lumA_clumped_female <- mr_pipeline_clumped_female(clumped_data = bmi_mrobject_women_clumped, 
                                                      outcome_data = lumA_data, 
                                                      riskfactor_label = "Body mass index (kg/m2)", 
                                                      outcome_label = "Luminal A-like")

bmi_lumBHER2neg_clumped_female <- mr_pipeline_clumped_female(clumped_data = bmi_mrobject_women_clumped, 
                                                             outcome_data = lumBHER2neg_data, 
                                                             riskfactor_label = "Body mass index (kg/m2)", 
                                                             outcome_label = "Luminal B/HER2-negative like")

bmi_lumB_clumped_female <- mr_pipeline_clumped_female(clumped_data = bmi_mrobject_women_clumped, 
                                                      outcome_data = lumB_data,
                                                      riskfactor_label = "Body mass index (kg/m2)", 
                                                      outcome_label = "Luminal B-like")

bmi_HER2enr_clumped_female <- mr_pipeline_clumped_female(clumped_data = bmi_mrobject_women_clumped, 
                                                         outcome_data = HER2enr_data, 
                                                         riskfactor_label = "Body mass index (kg/m2)", 
                                                         outcome_label = "HER2-enriched-like")

bmi_TNBC_clumped_female <- mr_pipeline_clumped_female(clumped_data = bmi_mrobject_women_clumped, 
                                                      outcome_data = TNBC_data, 
                                                      riskfactor_label = "Body mass index (kg/m2)", 
                                                      outcome_label = "Triple negative")

bmi_all_clumped_female <- rbind(bmi_overall_clumped_female,
                                bmi_lumA_clumped_female,
                                bmi_lumBHER2neg_clumped_female,
                                bmi_lumB_clumped_female,
                                bmi_HER2enr_clumped_female,
                                bmi_TNBC_clumped_female)

bmi_all_clumped_female

rm(bmi_overall_clumped_female, bmi_lumA_clumped_female, bmi_lumBHER2neg_clumped_female, bmi_lumB_clumped_female, bmi_HER2enr_clumped_female, bmi_TNBC_clumped_female)

rm(bmi_mrobject_women_clumped)


# Save bmi_all_clumped_female in separate Rdata object
# save(bmi_all_clumped_female, file = "manuscript/data_bmi_femalespecific_clumped.Rdata") # saved dd 01.12.2022 15:05

rm(bmi_all_clumped_female)



################################# T2D ##########################################

# Run pipeline for type 2 diabetes analyses including LD matrix

t2d_overall_ld <- mr_pipeline_ld(riskfactor_data = t2d_mrobject_def, 
                                 outcome_data = overallbc_data, 
                                 riskfactor_label = "Type 2 diabetes", 
                                 outcome_label = "All BCAC breast cancer cases",
                                 analysis = "IVW including LD matrix")

t2d_lumA_ld <- mr_pipeline_ld(riskfactor_data = t2d_mrobject_def, 
                              outcome_data = lumA_data, 
                              riskfactor_label = "Type 2 diabetes", 
                              outcome_label = "Luminal A-like",
                              analysis = "IVW including LD matrix")

t2d_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = t2d_mrobject_def, 
                                     outcome_data = lumBHER2neg_data, 
                                     riskfactor_label = "Type 2 diabetes", 
                                     outcome_label = "Luminal B/HER2-negative like",
                                     analysis = "IVW including LD matrix")

t2d_lumB_ld <- mr_pipeline_ld(riskfactor_data = t2d_mrobject_def, 
                              outcome_data = lumB_data, 
                              riskfactor_label = "Type 2 diabetes", 
                              outcome_label = "Luminal B-like",
                              analysis = "IVW including LD matrix")

t2d_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = t2d_mrobject_def, 
                                 outcome_data = HER2enr_data, 
                                 riskfactor_label = "Type 2 diabetes", 
                                 outcome_label = "HER2-enriched-like",
                                 analysis = "IVW including LD matrix")

t2d_TNBC_ld <- mr_pipeline_ld(riskfactor_data = t2d_mrobject_def, 
                              outcome_data = TNBC_data, 
                              riskfactor_label = "Type 2 diabetes", 
                              outcome_label = "Triple negative",
                              analysis = "IVW including LD matrix")


t2d_overall_pca <- mr_pipeline_pca(riskfactor_data = t2d_mrobject_def, 
                                   outcome_data = overallbc_data, 
                                   riskfactor_label = "Type 2 diabetes", 
                                   outcome_label = "All BCAC breast cancer cases",
                                   analysis = "IVW including PCs")

t2d_lumA_pca <- mr_pipeline_pca(riskfactor_data = t2d_mrobject_def, 
                                outcome_data = lumA_data, 
                                riskfactor_label = "Type 2 diabetes", 
                                outcome_label = "Luminal A-like",
                                analysis = "IVW including PCs")

t2d_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = t2d_mrobject_def, 
                                       outcome_data = lumBHER2neg_data, 
                                       riskfactor_label = "Type 2 diabetes", 
                                       outcome_label = "Luminal B/HER2-negative like",
                                       analysis = "IVW including PCs")

t2d_lumB_pca <- mr_pipeline_pca(riskfactor_data = t2d_mrobject_def, 
                                outcome_data = lumB_data, 
                                riskfactor_label = "Type 2 diabetes", 
                                outcome_label = "Luminal B-like",
                                analysis = "IVW including PCs")

t2d_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = t2d_mrobject_def, 
                                   outcome_data = HER2enr_data, 
                                   riskfactor_label = "Type 2 diabetes", 
                                   outcome_label = "HER2-enriched-like",
                                   analysis = "IVW including PCs")

t2d_TNBC_pca <- mr_pipeline_pca(riskfactor_data = t2d_mrobject_def, 
                                outcome_data = TNBC_data, 
                                riskfactor_label = "Type 2 diabetes", 
                                outcome_label = "Triple negative",
                                analysis = "IVW including PCs")


t2d_all_fig1 <- rbind(t2d_overall_ld,
                      t2d_overall_pca,
                      t2d_lumA_ld,
                      t2d_lumA_pca,
                      t2d_lumBHER2neg_ld,
                      t2d_lumBHER2neg_pca,
                      t2d_lumB_ld,
                      t2d_lumB_pca,
                      t2d_HER2enr_ld,
                      t2d_HER2enr_pca,
                      t2d_TNBC_ld,
                      t2d_TNBC_pca)

t2d_all_fig1

rm(t2d_overall_ld, t2d_lumA_ld, t2d_lumBHER2neg_ld, t2d_lumB_ld, t2d_HER2enr_ld, t2d_TNBC_ld,
   t2d_overall_pca, t2d_lumA_pca, t2d_lumBHER2neg_pca, t2d_lumB_pca, t2d_HER2enr_pca, t2d_TNBC_pca)

# Save t2d_all_fig1 in separate Rdata object
# cgwtools::resave(t2d_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")



# Run pipeline for type 2 diabetes with clumped variants
t2d_overall_clumped <- mr_pipeline_clumped(riskfactor_data = t2d_mrobject_def, 
                                           outcome_data = overallbc_data, 
                                           riskfactor_label = "Type 2 diabetes", 
                                           outcome_label = "All BCAC breast cancer cases")

t2d_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = t2d_mrobject_def, 
                                        outcome_data = lumA_data, 
                                        riskfactor_label = "Type 2 diabetes", 
                                        outcome_label = "Luminal A-like")

t2d_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = t2d_mrobject_def, 
                                               outcome_data = lumBHER2neg_data, 
                                               riskfactor_label = "Type 2 diabetes", 
                                               outcome_label = "Luminal B/HER2-negative like")

t2d_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = t2d_mrobject_def, 
                                        outcome_data = lumB_data,
                                        riskfactor_label = "Type 2 diabetes", 
                                        outcome_label = "Luminal B-like")

t2d_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = t2d_mrobject_def, 
                                           outcome_data = HER2enr_data, 
                                           riskfactor_label = "Type 2 diabetes", 
                                           outcome_label = "HER2-enriched-like")

t2d_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = t2d_mrobject_def, 
                                        outcome_data = TNBC_data, 
                                        riskfactor_label = "Type 2 diabetes", 
                                        outcome_label = "Triple negative")

t2d_all_table2 <- rbind(t2d_overall_clumped,
                        t2d_lumA_clumped,
                        t2d_lumBHER2neg_clumped,
                        t2d_lumB_clumped,
                        t2d_HER2enr_clumped,
                        t2d_TNBC_clumped)

t2d_all_table2

rm(t2d_overall_clumped, t2d_lumA_clumped, t2d_lumBHER2neg_clumped, t2d_lumB_clumped, t2d_HER2enr_clumped, t2d_TNBC_clumped)

rm(t2d_mrobject_def)



############################### age at menarche ################################

# Analyses for age at menarche

# Main analyses; including LD matrix
menarche_overall_ld <- mr_pipeline_ld(riskfactor_data = menarche_mrobject_def, 
                                      outcome_data = overallbc_data, 
                                      riskfactor_label = "Age at menarche (years)", 
                                      outcome_label = "All BCAC breast cancer cases",
                                      analysis = "IVW including LD matrix")

menarche_lumA_ld <- mr_pipeline_ld(riskfactor_data = menarche_mrobject_def, 
                                   outcome_data = lumA_data, 
                                   riskfactor_label = "Age at menarche (years)", 
                                   outcome_label = "Luminal A-like",
                                   analysis = "IVW including LD matrix")

menarche_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = menarche_mrobject_def, 
                                          outcome_data = lumBHER2neg_data, 
                                          riskfactor_label = "Age at menarche (years)", 
                                          outcome_label = "Luminal B/HER2-negative like",
                                          analysis = "IVW including LD matrix")

menarche_lumB_ld <- mr_pipeline_ld(riskfactor_data = menarche_mrobject_def, 
                                   outcome_data = lumB_data, 
                                   riskfactor_label = "Age at menarche (years)", 
                                   outcome_label = "Luminal B-like",
                                   analysis = "IVW including LD matrix")

menarche_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = menarche_mrobject_def, 
                                      outcome_data = HER2enr_data, 
                                      riskfactor_label = "Age at menarche (years)", 
                                      outcome_label = "HER2-enriched-like",
                                      analysis = "IVW including LD matrix")

menarche_TNBC_ld <- mr_pipeline_ld(riskfactor_data = menarche_mrobject_def, 
                                   outcome_data = TNBC_data, 
                                   riskfactor_label = "Age at menarche (years)", 
                                   outcome_label = "Triple negative",
                                   analysis = "IVW including LD matrix")


menarche_overall_pca <- mr_pipeline_pca(riskfactor_data = menarche_mrobject_def, 
                                        outcome_data = overallbc_data, 
                                        riskfactor_label = "Age at menarche (years)", 
                                        outcome_label = "All BCAC breast cancer cases",
                                        analysis = "IVW including PCs")

menarche_lumA_pca <- mr_pipeline_pca(riskfactor_data = menarche_mrobject_def, 
                                     outcome_data = lumA_data, 
                                     riskfactor_label = "Age at menarche (years)", 
                                     outcome_label = "Luminal A-like",
                                     analysis = "IVW including PCs")

menarche_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = menarche_mrobject_def, 
                                            outcome_data = lumBHER2neg_data, 
                                            riskfactor_label = "Age at menarche (years)", 
                                            outcome_label = "Luminal B/HER2-negative like",
                                            analysis = "IVW including PCs")

menarche_lumB_pca <- mr_pipeline_pca(riskfactor_data = menarche_mrobject_def, 
                                     outcome_data = lumB_data, 
                                     riskfactor_label = "Age at menarche (years)", 
                                     outcome_label = "Luminal B-like",
                                     analysis = "IVW including PCs")

menarche_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = menarche_mrobject_def, 
                                        outcome_data = HER2enr_data, 
                                        riskfactor_label = "Age at menarche (years)", 
                                        outcome_label = "HER2-enriched-like",
                                        analysis = "IVW including PCs")

menarche_TNBC_pca <- mr_pipeline_pca(riskfactor_data = menarche_mrobject_def, 
                                     outcome_data = TNBC_data, 
                                     riskfactor_label = "Age at menarche (years)", 
                                     outcome_label = "Triple negative",
                                     analysis = "IVW including PCs")


menarche_all_fig1 <- rbind(menarche_overall_ld,
                           menarche_overall_pca,
                           menarche_lumA_ld,
                           menarche_lumA_pca,
                           menarche_lumBHER2neg_ld,
                           menarche_lumBHER2neg_pca,
                           menarche_lumB_ld,
                           menarche_lumB_pca,
                           menarche_HER2enr_ld,
                           menarche_HER2enr_pca,
                           menarche_TNBC_ld,
                           menarche_TNBC_pca)

menarche_all_fig1

rm(menarche_overall_ld, menarche_lumA_ld, menarche_lumBHER2neg_ld, menarche_lumB_ld, menarche_HER2enr_ld, menarche_TNBC_ld,
   menarche_overall_pca, menarche_lumA_pca, menarche_lumBHER2neg_pca, menarche_lumB_pca, menarche_HER2enr_pca, menarche_TNBC_pca)

# Add menarche_all_fig1 to existing Rdata object
# cgwtools::resave(menarche_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Sensitivity analyses; clumped variants
menarche_overall_clumped <- mr_pipeline_clumped(riskfactor_data = menarche_mrobject_def, 
                                                outcome_data = overallbc_data, 
                                                riskfactor_label = "Age at menarche (years)", 
                                                outcome_label = "All BCAC breast cancer cases")

menarche_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = menarche_mrobject_def, 
                                             outcome_data = lumA_data, 
                                             riskfactor_label = "Age at menarche (years)", 
                                             outcome_label = "Luminal A-like")

menarche_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = menarche_mrobject_def, 
                                                    outcome_data = lumBHER2neg_data, 
                                                    riskfactor_label = "Age at menarche (years)", 
                                                    outcome_label = "Luminal B/HER2-negative like")

menarche_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = menarche_mrobject_def, 
                                             outcome_data = lumB_data,
                                             riskfactor_label = "Age at menarche (years)", 
                                             outcome_label = "Luminal B-like")

menarche_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = menarche_mrobject_def, 
                                                outcome_data = HER2enr_data, 
                                                riskfactor_label = "Age at menarche (years)", 
                                                outcome_label = "HER2-enriched-like")

menarche_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = menarche_mrobject_def, 
                                             outcome_data = TNBC_data, 
                                             riskfactor_label = "Age at menarche (years)", 
                                             outcome_label = "Triple negative")

menarche_all_table2 <- rbind(menarche_overall_clumped,
                             menarche_lumA_clumped,
                             menarche_lumBHER2neg_clumped,
                             menarche_lumB_clumped,
                             menarche_HER2enr_clumped,
                             menarche_TNBC_clumped)

menarche_all_table2

rm(menarche_overall_clumped, menarche_lumA_clumped, menarche_lumBHER2neg_clumped, menarche_lumB_clumped, menarche_HER2enr_clumped, menarche_TNBC_clumped)


# Multivariable MR for age at menarche - BMI (clumped variants)

# Clump variants from menarche_mrobject_def
clumped_menarche <- TwoSampleMR::clump_data(dat = menarche_mrobject_def,
                                            clump_kb = 1000,
                                            clump_r2 = 0.001,
                                            clump_p1 = 4.99e-08,
                                            clump_p2 = 4.99e-08,
                                            pop = "EUR")

# Search for this clumped dataset BMI sum stats for two data sets
# 1. Yengo et al (2018); ieu-b-40
BMI_yengo <- TwoSampleMR::extract_outcome_data(snps = clumped_menarche$SNP,
                                               outcomes = "ieu-b-40",
                                               proxies = T,
                                               rsq = 0.8,
                                               align_alleles = 1, # yes
                                               palindromes = 1, # yes
                                               maf_threshold = 0.3) # 216 SNPs

head(BMI_yengo)


# 2. UKB sumstats Ben Elsworth et al (2018); ukb-b-19953
BMI_UKB <- TwoSampleMR::extract_outcome_data(snps = clumped_menarche$SNP,
                                             outcomes = "ukb-b-19953",
                                             proxies = T,
                                             rsq = 0.8,
                                             align_alleles = 1, # yes
                                             palindromes = 1, # yes
                                             maf_threshold = 0.3) # 291 SNPs; 1 duplicate

head(BMI_UKB)


# Harmonize menarche and BMI data
# Yengo et al. BMI data
menarche_bmi_harmon_yengo <- TwoSampleMR::harmonise_data(exposure_dat = clumped_menarche,
                                                         outcome_dat = BMI_yengo,
                                                         action = 2) # 216 SNPs

table(menarche_bmi_harmon_yengo$mr_keep) # 8 SNPs should be removed before conducting analyses


menarche_bmi_harmon_yengo_def <- menarche_bmi_harmon_yengo %>%
  filter(mr_keep == "TRUE") # correct


# UKB BMI data
menarche_bmi_harmon_ukb <- TwoSampleMR::harmonise_data(exposure_dat = clumped_menarche,
                                                       outcome_dat = BMI_UKB,
                                                       action = 2) # 293 SNPs; 3 duplicates (rs235696)

table(menarche_bmi_harmon_ukb$mr_keep) # 281 would remain

# Subset to SNPs that should be kept + unique rsids 
menarche_bmi_harmon_ukb_def <- menarche_bmi_harmon_ukb %>%
  filter(mr_keep == "TRUE" & !duplicated(SNP)) # correct

# Check if there are still duplicated rsids in this set
# table(duplicated(menarche_bmi_harmon_ukb_def$SNP))

# Further investigate this SNP
# check <- menarche_bmi_harmon_ukb_def %>%
# filter(SNP == "rs235696")

# check # entries only differ on se.outcome and pvalue.outcome columns, but differences are very small; decision 08/08/2022 pick random one through additional filter in line 1049


# Function to perform multivariable MR analysis

mvmr <- function(bmi_data, outcome_dat, bmi_source, outcome_label){
  
  # Harmonize menarche and breast cancer outcome data
  menarche_outcome_harmon <- TwoSampleMR::harmonise_data(exposure_dat = clumped_menarche[which(clumped_menarche$SNP %in% bmi_data$SNP), ],
                                                         outcome_dat = outcome_dat,
                                                         action = 2) # checked; order of the variants is the same as in menarche_bmi_harmon!
  # effect alleles also appear to be the same!
  
  # Define variables to create multivariable MR input
  betaX <- as.matrix(bmi_data[c("beta.exposure", "beta.outcome")])
  
  betaY <- menarche_outcome_harmon$beta.outcome
  
  betaXse <- as.matrix(bmi_data[c("se.exposure", "se.outcome")])
  
  betaYse <- menarche_outcome_harmon$se.outcome
  
  snps <- menarche_outcome_harmon$SNP
  
  effect_allele <- menarche_outcome_harmon$effect_allele.outcome
  
  other_allele <- menarche_outcome_harmon$other_allele.outcome
  
  # Create multivariable MR input
  mvmr_data <- MendelianRandomization::mr_mvinput(bx = betaX,
                                                  bxse = betaXse,
                                                  by = betaY,
                                                  byse = betaYse,
                                                  exposure = c("Age at menarche", paste("BMI", bmi_source, sep = "_")),
                                                  outcome = outcome_label,
                                                  snps = snps,
                                                  effect_allele = effect_allele,
                                                  other_allele = other_allele)
  
  # Perform multivariable MR                     
  mvivw <- MendelianRandomization::mr_mvivw(object = mvmr_data,
                                            model = "random")
  
  
} # function is correct!


# Perform actual multivariable MR for each breast cancer outcome
# BMI from Yengo et al. --> also perform univariable MR for menarche with these 208 variants

mvivw_overall_yengo <- mvmr(bmi_data = menarche_bmi_harmon_yengo_def,
                            outcome_dat = overallbc_data,
                            bmi_source = "Yengo",
                            outcome_label = "Overall BC")

mvivw_lumA_yengo <- mvmr(bmi_data = menarche_bmi_harmon_yengo_def,
                         outcome_dat = lumA_data,
                         bmi_source = "Yengo",
                         outcome_label = "Luminal A-like")

mvivw_lumBHER2neg_yengo <- mvmr(bmi_data = menarche_bmi_harmon_yengo_def,
                                outcome_dat = lumBHER2neg_data,
                                bmi_source = "Yengo",
                                outcome_label = "Luminal B/HER2-negative like")

mvivw_lumB_yengo <- mvmr(bmi_data = menarche_bmi_harmon_yengo_def,
                         outcome_dat = lumB_data,
                         bmi_source = "Yengo",
                         outcome_label = "Luminal B-like")

mvivw_HER2enr_yengo <- mvmr(bmi_data = menarche_bmi_harmon_yengo_def,
                            outcome_dat = HER2enr_data,
                            bmi_source = "Yengo",
                            outcome_label = "HER2-enriched-like")

mvivw_TNBC_yengo <- mvmr(bmi_data = menarche_bmi_harmon_yengo_def,
                         outcome_dat = TNBC_data,
                         bmi_source = "Yengo",
                         outcome_label = "Triple negative")



# BMI from UKB --> also perform univariable MR for age at menarche with these 280 variants
mvivw_overall_ukb <- mvmr(bmi_data = menarche_bmi_harmon_ukb_def,
                          outcome_dat = overallbc_data,
                          bmi_source = "ukb",
                          outcome_label = "Overall BC")

mvivw_lumA_ukb <- mvmr(bmi_data = menarche_bmi_harmon_ukb_def,
                       outcome_dat = lumA_data,
                       bmi_source = "ukb",
                       outcome_label = "Luminal A-like")

mvivw_lumBHER2neg_ukb <- mvmr(bmi_data = menarche_bmi_harmon_ukb_def,
                              outcome_dat = lumBHER2neg_data,
                              bmi_source = "ukb",
                              outcome_label = "Luminal B/HER2-negative like")

mvivw_lumB_ukb <- mvmr(bmi_data = menarche_bmi_harmon_ukb_def,
                       outcome_dat = lumB_data,
                       bmi_source = "ukb",
                       outcome_label = "Luminal B-like")

mvivw_HER2enr_ukb <- mvmr(bmi_data = menarche_bmi_harmon_ukb_def,
                          outcome_dat = HER2enr_data,
                          bmi_source = "ukb",
                          outcome_label = "HER2-enriched-like")

mvivw_TNBC_ukb <- mvmr(bmi_data = menarche_bmi_harmon_ukb_def,
                       outcome_dat = TNBC_data,
                       bmi_source = "ukb",
                       outcome_label = "Triple negative")



# Load results clumped variants into R for making plot
load("manuscript/data_table2_dd27062022.Rdata")

# Select only estimates for age at menarche + only IVW estimates
menarche <- all_riskfactors_table2 %>% filter(`Risk factor` == "Age at menarche (years)")
menarche[8:23] <- NULL

# Make new variables to store multivariable estimates in
menarche$IVW.mvmr.yengo <- NA
menarche$se.mvmr.yengo <- NA
menarche$IVW.mvmr.ukb <- NA
menarche$se.mvmr.ukb <- NA

# Store multivarible MR estimates into lists + label results
# Results BMI Yengo et al.
mvivw_yengo <- list(mvivw_overall_yengo, 
                    mvivw_lumA_yengo,
                    mvivw_lumBHER2neg_yengo,
                    mvivw_lumB_yengo,
                    mvivw_HER2enr_yengo,
                    mvivw_TNBC_yengo)

names(mvivw_yengo) <- c("overallbc_data", "lumA_data", "lumBHER2neg_data", "lumB_data", "HER2enr_data", "TNBC_data")

# Results UKB
mvivw_ukb <- list(mvivw_overall_ukb, 
                  mvivw_lumA_ukb,
                  mvivw_lumBHER2neg_ukb,
                  mvivw_lumB_ukb,
                  mvivw_HER2enr_ukb,
                  mvivw_TNBC_ukb)

names(mvivw_ukb) <- c("overallbc_data", "lumA_data", "lumBHER2neg_data", "lumB_data", "HER2enr_data", "TNBC_data")

# Extract multivariable estimates and add to menarche dataframe
for (i in 1:6){
  
  menarche[i, "IVW.mvmr.yengo"] <- mvivw_yengo[[i]]@Estimate[1]
  menarche[i, "se.mvmr.yengo"] <- mvivw_yengo[[i]]@StdError[1]
  menarche[i, "IVW.mvmr.ukb"] <- mvivw_ukb[[i]]@Estimate[1]
  menarche[i, "se.mvmr.ukb"] <- mvivw_ukb[[i]]@StdError[1]
  
}

menarche # correct

# Calculate heterogeneity across multivariable MR estimates

subtypes <- menarche %>%
  filter(Outcome != "All BCAC breast cancer cases")

# Perform meta-analysis
meta_analysis.yengo <- rma(yi = subtypes$IVW.mvmr.yengo,
                           sei = subtypes$se.mvmr.yengo)

meta_analysis.yengo # 0% heterogeneity

# Perform meta-analysis
meta_analysis.ukb <- rma(yi = subtypes$IVW.mvmr.ukb,
                         sei = subtypes$se.mvmr.ukb)

meta_analysis.ukb # 0% heterogeneity


# Perform additional exclusion of estimate for triple negative breast cancer
subtypes_noTNBC <- subtypes %>%
  filter(Outcome != "Triple negative")

# Perform meta-analysis
meta_analysis_noTNBC.yengo <- rma(yi = subtypes_noTNBC$IVW.mvmr.yengo,
                                  sei = subtypes_noTNBC$se.mvmr.yengo)

meta_analysis_noTNBC.yengo # 15.08%

# Perform meta-analysis
meta_analysis_noTNBC.ukb <- rma(yi = subtypes_noTNBC$IVW.mvmr.ukb,
                                sei = subtypes_noTNBC$se.mvmr.ukb)

meta_analysis_noTNBC.ukb # 0%



################################# age at menopause #############################

# Analyses for age at menopause

# Main analyses; including LD matrix
menopause_overall_ld <- mr_pipeline_ld(riskfactor_data = menopause_mrobject_def, 
                                       outcome_data = overallbc_data, 
                                       riskfactor_label = "Age at menopause (years)", 
                                       outcome_label = "All BCAC breast cancer cases",
                                       analysis = "IVW including LD matrix")

menopause_lumA_ld <- mr_pipeline_ld(riskfactor_data = menopause_mrobject_def, 
                                    outcome_data = lumA_data, 
                                    riskfactor_label = "Age at menopause (years)", 
                                    outcome_label = "Luminal A-like",
                                    analysis = "IVW including LD matrix")

menopause_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = menopause_mrobject_def, 
                                           outcome_data = lumBHER2neg_data, 
                                           riskfactor_label = "Age at menopause (years)", 
                                           outcome_label = "Luminal B/HER2-negative like",
                                           analysis = "IVW including LD matrix")

menopause_lumB_ld <- mr_pipeline_ld(riskfactor_data = menopause_mrobject_def, 
                                    outcome_data = lumB_data, 
                                    riskfactor_label = "Age at menopause (years)", 
                                    outcome_label = "Luminal B-like",
                                    analysis = "IVW including LD matrix")

menopause_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = menopause_mrobject_def, 
                                       outcome_data = HER2enr_data, 
                                       riskfactor_label = "Age at menopause (years)", 
                                       outcome_label = "HER2-enriched-like",
                                       analysis = "IVW including LD matrix")

menopause_TNBC_ld <- mr_pipeline_ld(riskfactor_data = menopause_mrobject_def, 
                                    outcome_data = TNBC_data, 
                                    riskfactor_label = "Age at menopause (years)", 
                                    outcome_label = "Triple negative",
                                    analysis = "IVW including LD matrix")


menopause_overall_pca <- mr_pipeline_pca(riskfactor_data = menopause_mrobject_def, 
                                         outcome_data = overallbc_data, 
                                         riskfactor_label = "Age at menopause (years)", 
                                         outcome_label = "All BCAC breast cancer cases",
                                         analysis = "IVW including PCs")

menopause_lumA_pca <- mr_pipeline_pca(riskfactor_data = menopause_mrobject_def, 
                                      outcome_data = lumA_data, 
                                      riskfactor_label = "Age at menopause (years)", 
                                      outcome_label = "Luminal A-like",
                                      analysis = "IVW including PCs")

menopause_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = menopause_mrobject_def, 
                                             outcome_data = lumBHER2neg_data, 
                                             riskfactor_label = "Age at menopause (years)", 
                                             outcome_label = "Luminal B/HER2-negative like",
                                             analysis = "IVW including PCs")

menopause_lumB_pca <- mr_pipeline_pca(riskfactor_data = menopause_mrobject_def, 
                                      outcome_data = lumB_data, 
                                      riskfactor_label = "Age at menopause (years)", 
                                      outcome_label = "Luminal B-like",
                                      analysis = "IVW including PCs")

menopause_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = menopause_mrobject_def, 
                                         outcome_data = HER2enr_data, 
                                         riskfactor_label = "Age at menopause (years)", 
                                         outcome_label = "HER2-enriched-like",
                                         analysis = "IVW including PCs")

menopause_TNBC_pca <- mr_pipeline_pca(riskfactor_data = menopause_mrobject_def, 
                                      outcome_data = TNBC_data, 
                                      riskfactor_label = "Age at menopause (years)", 
                                      outcome_label = "Triple negative",
                                      analysis = "IVW including PCs")

menopause_all_fig1 <- rbind(menopause_overall_ld,
                            menopause_overall_pca,
                            menopause_lumA_ld,
                            menopause_lumA_pca,
                            menopause_lumBHER2neg_ld,
                            menopause_lumBHER2neg_pca,
                            menopause_lumB_ld,
                            menopause_lumB_pca,
                            menopause_HER2enr_ld,
                            menopause_HER2enr_pca,
                            menopause_TNBC_ld,
                            menopause_TNBC_pca)

menopause_all_fig1

rm(menopause_overall_ld, menopause_lumA_ld, menopause_lumBHER2neg_ld, menopause_lumB_ld, menopause_HER2enr_ld, menopause_TNBC_ld,
   menopause_overall_pca, menopause_lumA_pca, menopause_lumBHER2neg_pca, menopause_lumB_pca, menopause_HER2enr_pca, menopause_TNBC_pca)

# Add menopause_all_fig1 to existing Rdata object
# cgwtools::resave(menopause_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Sensitivity analyses; clumped variants

menopause_overall_clumped <- mr_pipeline_clumped(riskfactor_data = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], # otherwise the 2 NAs result in an error
                                                 outcome_data = overallbc_data, 
                                                 riskfactor_label = "Age at menopause (years)", 
                                                 outcome_label = "All BCAC breast cancer cases")

menopause_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], 
                                              outcome_data = lumA_data, 
                                              riskfactor_label = "Age at menopause (years)", 
                                              outcome_label = "Luminal A-like")

menopause_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], 
                                                     outcome_data = lumBHER2neg_data, 
                                                     riskfactor_label = "Age at menopause (years)", 
                                                     outcome_label = "Luminal B/HER2-negative like")

menopause_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], 
                                              outcome_data = lumB_data,
                                              riskfactor_label = "Age at menopause (years)", 
                                              outcome_label = "Luminal B-like")

menopause_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], 
                                                 outcome_data = HER2enr_data, 
                                                 riskfactor_label = "Age at menopause (years)", 
                                                 outcome_label = "HER2-enriched-like")

menopause_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], 
                                              outcome_data = TNBC_data, 
                                              riskfactor_label = "Age at menopause (years)", 
                                              outcome_label = "Triple negative")

menopause_all_table2 <- rbind(menopause_overall_clumped,
                              menopause_lumA_clumped,
                              menopause_lumBHER2neg_clumped,
                              menopause_lumB_clumped,
                              menopause_HER2enr_clumped,
                              menopause_TNBC_clumped)

menopause_all_table2

rm(menopause_overall_clumped, menopause_lumA_clumped, menopause_lumBHER2neg_clumped, menopause_lumB_clumped, menopause_HER2enr_clumped, menopause_TNBC_clumped)

rm(menopause_mrobject_def)



########################## percent breast density ##############################

# Analyses for percent breast density

# Main analyses; including LD matrix
density_overall_ld <- mr_pipeline_ld(riskfactor_data = density_mrobject_def, 
                                     outcome_data = overallbc_data, 
                                     riskfactor_label = "Breast density (%)", 
                                     outcome_label = "All BCAC breast cancer cases",
                                     analysis = "IVW including LD matrix")

density_lumA_ld <- mr_pipeline_ld(riskfactor_data = density_mrobject_def, 
                                  outcome_data = lumA_data, 
                                  riskfactor_label = "Breast density (%)", 
                                  outcome_label = "Luminal A-like",
                                  analysis = "IVW including LD matrix")

density_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = density_mrobject_def, 
                                         outcome_data = lumBHER2neg_data, 
                                         riskfactor_label = "Breast density (%)", 
                                         outcome_label = "Luminal B/HER2-negative like",
                                         analysis = "IVW including LD matrix")

density_lumB_ld <- mr_pipeline_ld(riskfactor_data = density_mrobject_def, 
                                  outcome_data = lumB_data, 
                                  riskfactor_label = "Breast density (%)", 
                                  outcome_label = "Luminal B-like",
                                  analysis = "IVW including LD matrix")

density_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = density_mrobject_def, 
                                     outcome_data = HER2enr_data, 
                                     riskfactor_label = "Breast density (%)", 
                                     outcome_label = "HER2-enriched-like",
                                     analysis = "IVW including LD matrix")

density_TNBC_ld <- mr_pipeline_ld(riskfactor_data = density_mrobject_def, 
                                  outcome_data = TNBC_data, 
                                  riskfactor_label = "Breast density (%)", 
                                  outcome_label = "Triple negative",
                                  analysis = "IVW including LD matrix")


density_overall_pca <- mr_pipeline_pca(riskfactor_data = density_mrobject_def, 
                                       outcome_data = overallbc_data, 
                                       riskfactor_label = "Breast density (%)", 
                                       outcome_label = "All BCAC breast cancer cases",
                                       analysis = "IVW including PCs")

density_lumA_pca <- mr_pipeline_pca(riskfactor_data = density_mrobject_def, 
                                    outcome_data = lumA_data, 
                                    riskfactor_label = "Breast density (%)", 
                                    outcome_label = "Luminal A-like",
                                    analysis = "IVW including PCs")

density_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = density_mrobject_def, 
                                           outcome_data = lumBHER2neg_data, 
                                           riskfactor_label = "Breast density (%)", 
                                           outcome_label = "Luminal B/HER2-negative like",
                                           analysis = "IVW including PCs")

density_lumB_pca <- mr_pipeline_pca(riskfactor_data = density_mrobject_def, 
                                    outcome_data = lumB_data, 
                                    riskfactor_label = "Breast density (%)", 
                                    outcome_label = "Luminal B-like",
                                    analysis = "IVW including PCs")

density_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = density_mrobject_def, 
                                       outcome_data = HER2enr_data, 
                                       riskfactor_label = "Breast density (%)", 
                                       outcome_label = "HER2-enriched-like",
                                       analysis = "IVW including PCs")

density_TNBC_pca <- mr_pipeline_pca(riskfactor_data = density_mrobject_def, 
                                    outcome_data = TNBC_data, 
                                    riskfactor_label = "Breast density (%)", 
                                    outcome_label = "Triple negative",
                                    analysis = "IVW including PCs")


density_all_fig1 <- rbind(density_overall_ld,
                          density_overall_pca,
                          density_lumA_ld,
                          density_lumA_pca,
                          density_lumBHER2neg_ld,
                          density_lumBHER2neg_pca,
                          density_lumB_ld,
                          density_lumB_pca,
                          density_HER2enr_ld,
                          density_HER2enr_pca,
                          density_TNBC_ld,
                          density_TNBC_pca)

density_all_fig1

rm(density_overall_ld, density_lumA_ld, density_lumBHER2neg_ld, density_lumB_ld, density_HER2enr_ld, density_TNBC_ld,
   density_overall_pca, density_lumA_pca, density_lumBHER2neg_pca, density_lumB_pca, density_HER2enr_pca, density_TNBC_pca)

# Add density_all_fig1 to existing Rdata object
# cgwtools::resave(density_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Sensitivity analyses; clumped variants
density_overall_clumped <- mr_pipeline_clumped(riskfactor_data = density_mrobject_def, 
                                               outcome_data = overallbc_data, 
                                               riskfactor_label = "Breast density (%)", 
                                               outcome_label = "All BCAC breast cancer cases") 

density_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = density_mrobject_def, 
                                            outcome_data = lumA_data, 
                                            riskfactor_label = "Breast density (%)", 
                                            outcome_label = "Luminal A-like")

density_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = density_mrobject_def, 
                                                   outcome_data = lumBHER2neg_data, 
                                                   riskfactor_label = "Breast density (%)", 
                                                   outcome_label = "Luminal B/HER2-negative like")

density_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = density_mrobject_def, 
                                            outcome_data = lumB_data,
                                            riskfactor_label = "Breast density (%)", 
                                            outcome_label = "Luminal B-like")

density_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = density_mrobject_def, 
                                               outcome_data = HER2enr_data, 
                                               riskfactor_label = "Breast density (%)", 
                                               outcome_label = "HER2-enriched-like")

density_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = density_mrobject_def, 
                                            outcome_data = TNBC_data, 
                                            riskfactor_label = "Breast density (%)", 
                                            outcome_label = "Triple negative")

density_all_table2 <- rbind(density_overall_clumped,
                            density_lumA_clumped,
                            density_lumBHER2neg_clumped,
                            density_lumB_clumped,
                            density_HER2enr_clumped,
                            density_TNBC_clumped)

density_all_table2

rm(density_overall_clumped, density_lumA_clumped, density_lumBHER2neg_clumped, density_lumB_clumped, density_HER2enr_clumped, density_TNBC_clumped)

rm(density_mrobject_def)



############################# alcohol consumption ##############################

# Analyses for alcohol consumption (drinks/week)

# Main analyses; including LD matrix
alcohol_overall_ld <- mr_pipeline_ld(riskfactor_data = alcohol_mrobject_def, 
                                     outcome_data = overallbc_data, 
                                     riskfactor_label = "Alcohol consumption (drinks/week)", 
                                     outcome_label = "All BCAC breast cancer cases",
                                     analysis = "IVW including LD matrix")

alcohol_lumA_ld <- mr_pipeline_ld(riskfactor_data = alcohol_mrobject_def, 
                                  outcome_data = lumA_data, 
                                  riskfactor_label = "Alcohol consumption (drinks/week)", 
                                  outcome_label = "Luminal A-like",
                                  analysis = "IVW including LD matrix")

alcohol_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = alcohol_mrobject_def, 
                                         outcome_data = lumBHER2neg_data, 
                                         riskfactor_label = "Alcohol consumption (drinks/week)", 
                                         outcome_label = "Luminal B/HER2-negative like",
                                         analysis = "IVW including LD matrix")

alcohol_lumB_ld <- mr_pipeline_ld(riskfactor_data = alcohol_mrobject_def, 
                                  outcome_data = lumB_data, 
                                  riskfactor_label = "Alcohol consumption (drinks/week)", 
                                  outcome_label = "Luminal B-like",
                                  analysis = "IVW including LD matrix")

alcohol_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = alcohol_mrobject_def, 
                                     outcome_data = HER2enr_data, 
                                     riskfactor_label = "Alcohol consumption (drinks/week)", 
                                     outcome_label = "HER2-enriched-like",
                                     analysis = "IVW including LD matrix")

alcohol_TNBC_ld <- mr_pipeline_ld(riskfactor_data = alcohol_mrobject_def, 
                                  outcome_data = TNBC_data, 
                                  riskfactor_label = "Alcohol consumption (drinks/week)", 
                                  outcome_label = "Triple negative",
                                  analysis = "IVW including LD matrix")


alcohol_overall_pca <- mr_pipeline_pca(riskfactor_data = alcohol_mrobject_def, 
                                       outcome_data = overallbc_data, 
                                       riskfactor_label = "Alcohol consumption (drinks/week)", 
                                       outcome_label = "All BCAC breast cancer cases",
                                       analysis = "IVW including PCs")

alcohol_lumA_pca <- mr_pipeline_pca(riskfactor_data = alcohol_mrobject_def, 
                                    outcome_data = lumA_data, 
                                    riskfactor_label = "Alcohol consumption (drinks/week)", 
                                    outcome_label = "Luminal A-like",
                                    analysis = "IVW including PCs")

alcohol_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = alcohol_mrobject_def, 
                                           outcome_data = lumBHER2neg_data, 
                                           riskfactor_label = "Alcohol consumption (drinks/week)", 
                                           outcome_label = "Luminal B/HER2-negative like",
                                           analysis = "IVW including PCs")

alcohol_lumB_pca <- mr_pipeline_pca(riskfactor_data = alcohol_mrobject_def, 
                                    outcome_data = lumB_data, 
                                    riskfactor_label = "Alcohol consumption (drinks/week)", 
                                    outcome_label = "Luminal B-like",
                                    analysis = "IVW including PCs")

alcohol_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = alcohol_mrobject_def, 
                                       outcome_data = HER2enr_data, 
                                       riskfactor_label = "Alcohol consumption (drinks/week)", 
                                       outcome_label = "HER2-enriched-like",
                                       analysis = "IVW including PCs")

alcohol_TNBC_pca <- mr_pipeline_pca(riskfactor_data = alcohol_mrobject_def, 
                                    outcome_data = TNBC_data, 
                                    riskfactor_label = "Alcohol consumption (drinks/week)", 
                                    outcome_label = "Triple negative",
                                    analysis = "IVW including PCs")


alcohol_all_fig1 <- rbind(alcohol_overall_ld,
                          alcohol_overall_pca,
                          alcohol_lumA_ld,
                          alcohol_lumA_pca,
                          alcohol_lumBHER2neg_ld,
                          alcohol_lumBHER2neg_pca,
                          alcohol_lumB_ld,
                          alcohol_lumB_pca,
                          alcohol_HER2enr_ld,
                          alcohol_HER2enr_pca,
                          alcohol_TNBC_ld,
                          alcohol_TNBC_pca)

alcohol_all_fig1

rm(alcohol_overall_ld, alcohol_lumA_ld, alcohol_lumBHER2neg_ld, alcohol_lumB_ld, alcohol_HER2enr_ld, alcohol_TNBC_ld,
   alcohol_overall_pca, alcohol_lumA_pca, alcohol_lumBHER2neg_pca, alcohol_lumB_pca, alcohol_HER2enr_pca, alcohol_TNBC_pca)


# Add alcohol_all_fig1 to existing Rdata object
# cgwtools::resave(alcohol_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Sensitivity analyses; clumped variants
alcohol_overall_clumped <- mr_pipeline_clumped(riskfactor_data = alcohol_mrobject_def, 
                                               outcome_data = overallbc_data, 
                                               riskfactor_label = "Alcohol consumption (drinks/week)", 
                                               outcome_label = "All BCAC breast cancer cases")

alcohol_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = alcohol_mrobject_def, 
                                            outcome_data = lumA_data, 
                                            riskfactor_label = "Alcohol consumption (drinks/week)", 
                                            outcome_label = "Luminal A-like")

alcohol_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = alcohol_mrobject_def, 
                                                   outcome_data = lumBHER2neg_data, 
                                                   riskfactor_label = "Alcohol consumption (drinks/week)", 
                                                   outcome_label = "Luminal B/HER2-negative like")

alcohol_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = alcohol_mrobject_def, 
                                            outcome_data = lumB_data,
                                            riskfactor_label = "Alcohol consumption (drinks/week)", 
                                            outcome_label = "Luminal B-like")

alcohol_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = alcohol_mrobject_def, 
                                               outcome_data = HER2enr_data, 
                                               riskfactor_label = "Alcohol consumption (drinks/week)", 
                                               outcome_label = "HER2-enriched-like")

alcohol_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = alcohol_mrobject_def, 
                                            outcome_data = TNBC_data, 
                                            riskfactor_label = "Alcohol consumption (drinks/week)", 
                                            outcome_label = "Triple negative")

alcohol_all_table2 <- rbind(alcohol_overall_clumped,
                            alcohol_lumA_clumped,
                            alcohol_lumBHER2neg_clumped,
                            alcohol_lumB_clumped,
                            alcohol_HER2enr_clumped,
                            alcohol_TNBC_clumped)

alcohol_all_table2

rm(alcohol_overall_clumped, alcohol_lumA_clumped, alcohol_lumBHER2neg_clumped, alcohol_lumB_clumped, alcohol_HER2enr_clumped, alcohol_TNBC_clumped)

rm(alcohol_mrobject_def)



################################ smoking #######################################

# Analyses for smoking behaviour (ever smoked regularly; yes/no)

# Main analyses; including LD matrix
smoking_overall_ld <- mr_pipeline_ld(riskfactor_data = smoking_mrobject_def, 
                                     outcome_data = overallbc_data, 
                                     riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                     outcome_label = "All BCAC breast cancer cases",
                                     analysis = "IVW including LD matrix")

smoking_lumA_ld <- mr_pipeline_ld(riskfactor_data = smoking_mrobject_def, 
                                  outcome_data = lumA_data, 
                                  riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                  outcome_label = "Luminal A-like",
                                  analysis = "IVW including LD matrix")

smoking_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = smoking_mrobject_def, 
                                         outcome_data = lumBHER2neg_data, 
                                         riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                         outcome_label = "Luminal B/HER2-negative like",
                                         analysis = "IVW including LD matrix")

smoking_lumB_ld <- mr_pipeline_ld(riskfactor_data = smoking_mrobject_def, 
                                  outcome_data = lumB_data, 
                                  riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                  outcome_label = "Luminal B-like",
                                  analysis = "IVW including LD matrix")

smoking_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = smoking_mrobject_def, 
                                     outcome_data = HER2enr_data, 
                                     riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                     outcome_label = "HER2-enriched-like",
                                     analysis = "IVW including LD matrix")

smoking_TNBC_ld <- mr_pipeline_ld(riskfactor_data = smoking_mrobject_def, 
                                  outcome_data = TNBC_data, 
                                  riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                  outcome_label = "Triple negative",
                                  analysis = "IVW including LD matrix")


smoking_overall_pca <- mr_pipeline_pca(riskfactor_data = smoking_mrobject_def, 
                                       outcome_data = overallbc_data, 
                                       riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                       outcome_label = "All BCAC breast cancer cases",
                                       analysis = "IVW including PCs")

smoking_lumA_pca <- mr_pipeline_pca(riskfactor_data = smoking_mrobject_def, 
                                    outcome_data = lumA_data, 
                                    riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                    outcome_label = "Luminal A-like",
                                    analysis = "IVW including PCs")

smoking_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = smoking_mrobject_def, 
                                           outcome_data = lumBHER2neg_data, 
                                           riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                           outcome_label = "Luminal B/HER2-negative like",
                                           analysis = "IVW including PCs")

smoking_lumB_pca <- mr_pipeline_pca(riskfactor_data = smoking_mrobject_def, 
                                    outcome_data = lumB_data, 
                                    riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                    outcome_label = "Luminal B-like",
                                    analysis = "IVW including PCs")

smoking_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = smoking_mrobject_def, 
                                       outcome_data = HER2enr_data, 
                                       riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                       outcome_label = "HER2-enriched-like",
                                       analysis = "IVW including PCs")

smoking_TNBC_pca <- mr_pipeline_pca(riskfactor_data = smoking_mrobject_def, 
                                    outcome_data = TNBC_data, 
                                    riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                    outcome_label = "Triple negative",
                                    analysis = "IVW including PCs")


smoking_all_fig1 <- rbind(smoking_overall_ld,
                          smoking_overall_pca,
                          smoking_lumA_ld,
                          smoking_lumA_pca,
                          smoking_lumBHER2neg_ld,
                          smoking_lumBHER2neg_pca,
                          smoking_lumB_ld,
                          smoking_lumB_pca,
                          smoking_HER2enr_ld,
                          smoking_HER2enr_pca,
                          smoking_TNBC_ld,
                          smoking_TNBC_pca)

smoking_all_fig1

rm(smoking_overall_ld, smoking_lumA_ld, smoking_lumBHER2neg_ld, smoking_lumB_ld, smoking_HER2enr_ld, smoking_TNBC_ld,
   smoking_overall_pca, smoking_lumA_pca, smoking_lumBHER2neg_pca, smoking_lumB_pca, smoking_HER2enr_pca, smoking_TNBC_pca)

# Add smoking_all_fig1 to existing Rdata object
# cgwtools::resave(smoking_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Sensitivity analyses; clumped variants
smoking_overall_clumped <- mr_pipeline_clumped(riskfactor_data = smoking_mrobject_def, 
                                               outcome_data = overallbc_data, 
                                               riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                               outcome_label = "All BCAC breast cancer cases")

smoking_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = smoking_mrobject_def, 
                                            outcome_data = lumA_data, 
                                            riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                            outcome_label = "Luminal A-like")

smoking_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = smoking_mrobject_def, 
                                                   outcome_data = lumBHER2neg_data, 
                                                   riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                                   outcome_label = "Luminal B/HER2-negative like")

smoking_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = smoking_mrobject_def, 
                                            outcome_data = lumB_data,
                                            riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                            outcome_label = "Luminal B-like")

smoking_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = smoking_mrobject_def, 
                                               outcome_data = HER2enr_data, 
                                               riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                               outcome_label = "HER2-enriched-like")

smoking_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = smoking_mrobject_def, 
                                            outcome_data = TNBC_data, 
                                            riskfactor_label = "Smoking (ever smoked regularly; yes/no)", 
                                            outcome_label = "Triple negative")

smoking_all_table2 <- rbind(smoking_overall_clumped,
                            smoking_lumA_clumped,
                            smoking_lumBHER2neg_clumped,
                            smoking_lumB_clumped,
                            smoking_HER2enr_clumped,
                            smoking_TNBC_clumped)

smoking_all_table2

rm(smoking_overall_clumped, smoking_lumA_clumped, smoking_lumBHER2neg_clumped, smoking_lumB_clumped, smoking_HER2enr_clumped, smoking_TNBC_clumped)

rm(smoking_mrobject_def)



################################## physical activity ###########################

# Analyses for physical activity

# Main analyses; including LD matrix
activity_overall_ld <- mr_pipeline_ld(riskfactor_data = activity_mrobject_def, 
                                      outcome_data = overallbc_data, 
                                      riskfactor_label = "Overall activity", 
                                      outcome_label = "All BCAC breast cancer cases",
                                      analysis = "IVW including LD matrix")

activity_lumA_ld <- mr_pipeline_ld(riskfactor_data = activity_mrobject_def, 
                                   outcome_data = lumA_data, 
                                   riskfactor_label = "Overall activity", 
                                   outcome_label = "Luminal A-like",
                                   analysis = "IVW including LD matrix")

activity_lumBHER2neg_ld <- mr_pipeline_ld(riskfactor_data = activity_mrobject_def, 
                                          outcome_data = lumBHER2neg_data, 
                                          riskfactor_label = "Overall activity", 
                                          outcome_label = "Luminal B/HER2-negative like",
                                          analysis = "IVW including LD matrix")

activity_lumB_ld <- mr_pipeline_ld(riskfactor_data = activity_mrobject_def, 
                                   outcome_data = lumB_data, 
                                   riskfactor_label = "Overall activity", 
                                   outcome_label = "Luminal B-like",
                                   analysis = "IVW including LD matrix")

activity_HER2enr_ld <- mr_pipeline_ld(riskfactor_data = activity_mrobject_def, 
                                      outcome_data = HER2enr_data, 
                                      riskfactor_label = "Overall activity", 
                                      outcome_label = "HER2-enriched-like",
                                      analysis = "IVW including LD matrix")

activity_TNBC_ld <- mr_pipeline_ld(riskfactor_data = activity_mrobject_def, 
                                   outcome_data = TNBC_data, 
                                   riskfactor_label = "Overall activity", 
                                   outcome_label = "Triple negative",
                                   analysis = "IVW including LD matrix")


activity_overall_pca <- mr_pipeline_pca(riskfactor_data = activity_mrobject_def, 
                                        outcome_data = overallbc_data, 
                                        riskfactor_label = "Overall activity", 
                                        outcome_label = "All BCAC breast cancer cases",
                                        analysis = "IVW including PCs")

activity_lumA_pca <- mr_pipeline_pca(riskfactor_data = activity_mrobject_def, 
                                     outcome_data = lumA_data, 
                                     riskfactor_label = "Overall activity", 
                                     outcome_label = "Luminal A-like",
                                     analysis = "IVW including PCs")

activity_lumBHER2neg_pca <- mr_pipeline_pca(riskfactor_data = activity_mrobject_def, 
                                            outcome_data = lumBHER2neg_data, 
                                            riskfactor_label = "Overall activity", 
                                            outcome_label = "Luminal B/HER2-negative like",
                                            analysis = "IVW including PCs")

activity_lumB_pca <- mr_pipeline_pca(riskfactor_data = activity_mrobject_def, 
                                     outcome_data = lumB_data, 
                                     riskfactor_label = "Overall activity", 
                                     outcome_label = "Luminal B-like",
                                     analysis = "IVW including PCs")

activity_HER2enr_pca <- mr_pipeline_pca(riskfactor_data = activity_mrobject_def, 
                                        outcome_data = HER2enr_data, 
                                        riskfactor_label = "Overall activity", 
                                        outcome_label = "HER2-enriched-like",
                                        analysis = "IVW including PCs")

activity_TNBC_pca <- mr_pipeline_pca(riskfactor_data = activity_mrobject_def, 
                                     outcome_data = TNBC_data, 
                                     riskfactor_label = "Overall activity", 
                                     outcome_label = "Triple negative",
                                     analysis = "IVW including PCs")


activity_all_fig1 <- rbind(activity_overall_ld,
                           activity_overall_pca,
                           activity_lumA_ld,
                           activity_lumA_pca,
                           activity_lumBHER2neg_ld,
                           activity_lumBHER2neg_pca,
                           activity_lumB_ld,
                           activity_lumB_pca,
                           activity_HER2enr_ld,
                           activity_HER2enr_pca,
                           activity_TNBC_ld,
                           activity_TNBC_pca)

activity_all_fig1

rm(activity_overall_ld, activity_lumA_ld, activity_lumBHER2neg_ld, activity_lumB_ld, activity_HER2enr_ld, activity_TNBC_ld,
   activity_overall_pca, activity_lumA_pca, activity_lumBHER2neg_pca, activity_lumB_pca, activity_HER2enr_pca, activity_TNBC_pca)

# Add activity_all_fig1 to existing Rdata object
# cgwtools::resave(activity_all_fig1, file = "manuscript/data_fig1_dd02022022.Rdata")


# Sensitivity analyses; clumped variants
activity_overall_clumped <- mr_pipeline_clumped(riskfactor_data = activity_mrobject_def, 
                                                outcome_data = overallbc_data, 
                                                riskfactor_label = "Overall activity", 
                                                outcome_label = "All BCAC breast cancer cases")

activity_lumA_clumped <- mr_pipeline_clumped(riskfactor_data = activity_mrobject_def, 
                                             outcome_data = lumA_data, 
                                             riskfactor_label = "Overall activity", 
                                             outcome_label = "Luminal A-like")

activity_lumBHER2neg_clumped <- mr_pipeline_clumped(riskfactor_data = activity_mrobject_def, 
                                                    outcome_data = lumBHER2neg_data, 
                                                    riskfactor_label = "Overall activity", 
                                                    outcome_label = "Luminal B/HER2-negative like")

activity_lumB_clumped <- mr_pipeline_clumped(riskfactor_data = activity_mrobject_def, 
                                             outcome_data = lumB_data,
                                             riskfactor_label = "Overall activity", 
                                             outcome_label = "Luminal B-like")

activity_HER2enr_clumped <- mr_pipeline_clumped(riskfactor_data = activity_mrobject_def, 
                                                outcome_data = HER2enr_data, 
                                                riskfactor_label = "Overall activity", 
                                                outcome_label = "HER2-enriched-like")

activity_TNBC_clumped <- mr_pipeline_clumped(riskfactor_data = activity_mrobject_def, 
                                             outcome_data = TNBC_data, 
                                             riskfactor_label = "Overall activity", 
                                             outcome_label = "Triple negative")

activity_all_table2 <- rbind(activity_overall_clumped,
                             activity_lumA_clumped,
                             activity_lumBHER2neg_clumped,
                             activity_lumB_clumped,
                             activity_HER2enr_clumped,
                             activity_TNBC_clumped)

activity_all_table2

rm(activity_overall_clumped, activity_lumA_clumped, activity_lumBHER2neg_clumped, activity_lumB_clumped, activity_HER2enr_clumped, activity_TNBC_clumped)


# Run primary analyses for activity using female-specific weights for genetic instruments

# Load female specific MR object into R
load("sumstat_riskfactors/female_betas/activity_mrobject_women.RData")


activity_overall_pca_female <- mr_pipeline_pca(riskfactor_data = activity_mrobject_women, 
                                               outcome_data = overallbc_data, 
                                               riskfactor_label = "Overall activity", 
                                               outcome_label = "All BCAC breast cancer cases",
                                               analysis = "IVW including PCs; female betas")

activity_lumA_pca_female <- mr_pipeline_pca(riskfactor_data = activity_mrobject_women, 
                                            outcome_data = lumA_data, 
                                            riskfactor_label = "Overall activity", 
                                            outcome_label = "Luminal A-like",
                                            analysis = "IVW including PCs; female betas")

activity_lumBHER2neg_pca_female <- mr_pipeline_pca(riskfactor_data = activity_mrobject_women, 
                                                   outcome_data = lumBHER2neg_data, 
                                                   riskfactor_label = "Overall activity", 
                                                   outcome_label = "Luminal B/HER2-negative like",
                                                   analysis = "IVW including PCs; female betas")

activity_lumB_pca_female <- mr_pipeline_pca(riskfactor_data = activity_mrobject_women, 
                                            outcome_data = lumB_data, 
                                            riskfactor_label = "Overall activity", 
                                            outcome_label = "Luminal B-like",
                                            analysis = "IVW including PCs; female betas")

activity_HER2enr_pca_female <- mr_pipeline_pca(riskfactor_data = activity_mrobject_women, 
                                               outcome_data = HER2enr_data, 
                                               riskfactor_label = "Overall activity", 
                                               outcome_label = "HER2-enriched-like",
                                               analysis = "IVW including PCs; female betas")

activity_TNBC_pca_female <- mr_pipeline_pca(riskfactor_data = activity_mrobject_women, 
                                            outcome_data = TNBC_data, 
                                            riskfactor_label = "Overall activity", 
                                            outcome_label = "Triple negative",
                                            analysis = "IVW including PCs; female betas")


activity_pca_all_female <- rbind(activity_overall_pca_female,
                                 activity_lumA_pca_female,
                                 activity_lumBHER2neg_pca_female,
                                 activity_lumB_pca_female,
                                 activity_HER2enr_pca_female,
                                 activity_TNBC_pca_female)

activity_pca_all_female

rm(activity_overall_pca_female, activity_lumA_pca_female, activity_lumBHER2neg_pca_female, activity_lumB_pca_female, activity_HER2enr_pca_female, activity_TNBC_pca_female)


# Save activity_pca_all_female in separate Rdata object
# save(activity_pca_all_female, file = "manuscript/data_activity_femalespecific.Rdata") # saved 01-12-2022 at 15:11

rm(activity_pca_all_female)


# Repeat analyses with clumped variants with female-specific estimates
# Maybe best to do clumping based on sex-combined data; use that dataframe to subset female-specific dataframe
# and subsequently run the pipeline without clumping; because that has then already been "done"

activity_mrobject_clumped <- TwoSampleMR::clump_data(dat = activity_mrobject_def,
                                                     clump_kb = 1000,
                                                     clump_r2 = 0.001,
                                                     clump_p1 = 4.99e-08,
                                                     clump_p2 = 4.99e-08,
                                                     pop = "EUR")


activity_mrobject_women_clumped <- activity_mrobject_women %>%
  filter(SNP %in% activity_mrobject_clumped$SNP)

rm(activity_mrobject_clumped,
   activity_mrobject_def,
   activity_mrobject_women)


activity_overall_clumped_female <- mr_pipeline_clumped_female(clumped_data = activity_mrobject_women_clumped, 
                                                              outcome_data = overallbc_data, 
                                                              riskfactor_label = "Overall activity", 
                                                              outcome_label = "All BCAC breast cancer cases")

activity_lumA_clumped_female <- mr_pipeline_clumped_female(clumped_data = activity_mrobject_women_clumped, 
                                                           outcome_data = lumA_data, 
                                                           riskfactor_label = "Overall activity", 
                                                           outcome_label = "Luminal A-like")

activity_lumBHER2neg_clumped_female <- mr_pipeline_clumped_female(clumped_data = activity_mrobject_women_clumped, 
                                                                  outcome_data = lumBHER2neg_data, 
                                                                  riskfactor_label = "Overall activity", 
                                                                  outcome_label = "Luminal B/HER2-negative like")

activity_lumB_clumped_female <- mr_pipeline_clumped_female(clumped_data = activity_mrobject_women_clumped, 
                                                           outcome_data = lumB_data,
                                                           riskfactor_label = "Overall activity", 
                                                           outcome_label = "Luminal B-like")

activity_HER2enr_clumped_female <- mr_pipeline_clumped_female(clumped_data = activity_mrobject_women_clumped, 
                                                              outcome_data = HER2enr_data, 
                                                              riskfactor_label = "Overall activity", 
                                                              outcome_label = "HER2-enriched-like")

activity_TNBC_clumped_female <- mr_pipeline_clumped_female(clumped_data = activity_mrobject_women_clumped, 
                                                           outcome_data = TNBC_data, 
                                                           riskfactor_label = "Overall activity", 
                                                           outcome_label = "Triple negative")

activity_all_clumped_female <- rbind(activity_overall_clumped_female,
                                     activity_lumA_clumped_female,
                                     activity_lumBHER2neg_clumped_female,
                                     activity_lumB_clumped_female,
                                     activity_HER2enr_clumped_female,
                                     activity_TNBC_clumped_female)

activity_all_clumped_female

rm(activity_overall_clumped_female, activity_lumA_clumped_female, activity_lumBHER2neg_clumped_female, activity_lumB_clumped_female, activity_HER2enr_clumped_female, activity_TNBC_clumped_female)

rm(activity_mrobject_women_clumped)


# Save activity_all_clumped_female in separate Rdata object
# save(activity_all_clumped_female, file = "manuscript/data_activity_femalespecific_clumped.Rdata") # saved dd 01.12.2022 15:12

rm(activity_all_clumped_female)



################################################################################

# Calculate % of heterogeneity across subtype-specific estimates of the primary analysis

heterogeneity <- function(data){
  
  # Create dataframe to store results in
  I2_results = as.data.frame(matrix(nrow = 1, ncol = 3))
  
  colnames(I2_results) <- c("Risk factor", "I.squared", "I.squared.noTNBC")
  
  # Subset data to only estimates of IVW analyses including PCs and exclude estimate for overall breast cancer
  subtypes <- data %>%
    filter(Analysis == "IVW including PCs" & Outcome != "All BCAC breast cancer cases")
  
  # Perform meta-analysis
  meta_analysis <- rma(yi = subtypes$IVW.est,
                       sei = subtypes$se.est)
  
  # Save results in I2_results
  I2_results[1,1] = unique(subtypes$`Risk factor`)
  
  I2_results[1,2] = round(meta_analysis$I2, digits = 1)
  
  
  # Perform additional exlcusion of estimate for triple negative breast cancer
  subtypes_noTNBC <- subtypes %>%
    filter(Outcome != "Triple negative")
  
  # Perform meta-analysis in this leave-TNBC-out analysis
  meta_analysis_noTNBC <- rma(yi = subtypes_noTNBC$IVW.est,
                              sei = subtypes_noTNBC$se.est)
  
  # Store results also in I2_results
  I2_results[1,3] = round(meta_analysis_noTNBC$I2, digits = 1)
  
  
  return(I2_results)
  
}

# Height
heterogeneity_height <- heterogeneity(data = height_all_fig1)

# BMI
heterogeneity_bmi <- heterogeneity(data = bmi_all_fig1)

# T2D
heterogeneity_t2d <- heterogeneity(data = t2d_all_fig1)

# Age at menarche
heterogeneity_menarche <- heterogeneity(data = menarche_all_fig1)

# Age at menopause
heterogeneity_menopause <- heterogeneity(data = menopause_all_fig1)

# Percent density
heterogeneity_density <- heterogeneity(data = density_all_fig1)

# Alcohol consumption
heterogeneity_alcohol <- heterogeneity(data = alcohol_all_fig1)

# Smoking
heterogeneity_smoking <- heterogeneity(data = smoking_all_fig1)

# Physical activity
heterogeneity_activity <- heterogeneity(data = activity_all_fig1)

# Combine all estimates
all_I2 <- rbind(heterogeneity_height,
                heterogeneity_bmi,
                heterogeneity_t2d,
                heterogeneity_menarche,
                heterogeneity_menopause,
                heterogeneity_density,
                heterogeneity_alcohol,
                heterogeneity_smoking,
                heterogeneity_activity)

all_I2

# Remove individual I2 estimates
rm(heterogeneity_activity, heterogeneity_alcohol, heterogeneity_bmi, heterogeneity_density, heterogeneity_height,
   heterogeneity_menarche, heterogeneity_menopause, heterogeneity_smoking, heterogeneity_t2d)



################################################################################

# Calculate % of heterogeneity across subtype-specific estimates of the secondary analysis

# First make new MR presso estimate and se variable that includes IVW estimate and se from clumped analysis if MR-PRESSO est and se are NA
# The reason for this is that for these instances, no outliers were detected and thus the IVW estimate and se are assumed to be unbiased

heterogeneity_robust <- all_riskfactors_table2_inclpresso %>%
  mutate(PRESSO.est = case_when(
    !is.na(PRESSO.est.y) ~ PRESSO.est.y,
    is.na(PRESSO.est.y) ~ IVW.est)) %>%
  mutate(PRESSO.se = case_when(
    !is.na(PRESSO.se.y) ~ PRESSO.se.y,
    is.na(PRESSO.se.y) ~ IVW.se)) # checked = correct

# Adapt function to do same calculations for IVW analysis restricted to clumped variants + MR-Egger + Weighted median + weighted mode (dd 27.06.2022)

heterogeneity <- function(risk_factor){
  
  # Create dataframe to store results in
  I2_results = as.data.frame(matrix(nrow = 1, ncol = 11))
  
  colnames(I2_results) <- c("Risk factor", "I.squared.ivwclump", "I.squared.ivw.clump.noTNBC",
                            "I.squared.egger", "I.squared.egger.noTNBC",
                            "I.squared.median", "I.squared.median.noTNBC",
                            "I.squared.mode", "I.squared.mode.noTNBC",
                            "I.squared.presso", "I.squared.presso.noTNBC")
  
  # Subset data to only estimates of IVW analyses including PCs and exclude estimate for overall breast cancer
  subtypes <- heterogeneity_robust %>%
    filter(`Risk factor` == risk_factor & Outcome != "All BCAC breast cancer cases")
  
  # Perform meta-analysis
  # IVW
  meta_analysis.ivw <- rma(yi = subtypes$IVW.est,
                           sei = subtypes$IVW.se)
  
  # Save results in I2_results
  I2_results[1,1] = unique(subtypes$`Risk factor`)
  
  I2_results[1,2] = round(meta_analysis.ivw$I2, digits = 1)
  
  
  # MR-Egger
  meta_analysis.egger <- rma(yi = subtypes$Egger.est,
                             sei = subtypes$Egger.se)
  
  # Save results in I2_results
  I2_results[1,4] = round(meta_analysis.egger$I2, digits = 1)
  
  
  # Weighted median
  meta_analysis.median <- rma(yi = subtypes$Median.est,
                              sei = subtypes$Median.se)
  
  # Save results in I2_results
  I2_results[1,6] = round(meta_analysis.median$I2, digits = 1)
  
  
  # Weighted mode
  meta_analysis.mode <- rma(yi = subtypes$Mode.est,
                            sei = subtypes$Mode.se)
  
  # Save results in I2_results
  I2_results[1,8] = round(meta_analysis.mode$I2, digits = 1)
  
  
  # MR-PRESSO
  meta_analysis.presso <- rma(yi = subtypes$PRESSO.est,
                              sei = subtypes$PRESSO.se)
  
  # Save results in I2_results
  I2_results[1,10] = round(meta_analysis.presso$I2, digits = 1)
  
  
  
  # Perform additional exlcusion of estimate for triple negative breast cancer for each method
  subtypes_noTNBC <- subtypes %>%
    filter(Outcome != "Triple negative")
  
  # Perform meta-analysis in this leave-TNBC-out analysis
  # IVW
  meta_analysis_noTNBC.ivw <- rma(yi = subtypes_noTNBC$IVW.est,
                                  sei = subtypes_noTNBC$IVW.se)
  
  # Store results also in I2_results
  I2_results[1,3] = round(meta_analysis_noTNBC.ivw$I2, digits = 1)
  
  # MR-Egger
  meta_analysis_noTNBC.egger <- rma(yi = subtypes_noTNBC$Egger.est,
                                    sei = subtypes_noTNBC$Egger.se)
  
  # Store results also in I2_results
  I2_results[1,5] = round(meta_analysis_noTNBC.egger$I2, digits = 1)
  
  
  # Weighted median
  meta_analysis_noTNBC.median <- rma(yi = subtypes_noTNBC$Median.est,
                                     sei = subtypes_noTNBC$Median.se)
  
  # Store results also in I2_results
  I2_results[1,7] = round(meta_analysis_noTNBC.median$I2, digits = 1)
  
  
  # Weighted mode
  meta_analysis_noTNBC.mode <- rma(yi = subtypes_noTNBC$Mode.est,
                                   sei = subtypes_noTNBC$Mode.se)
  
  # Store results also in I2_results
  I2_results[1,9] = round(meta_analysis_noTNBC.mode$I2, digits = 1)
  
  
  # MR-PRESSO
  meta_analysis_noTNBC.presso <- rma(yi = subtypes_noTNBC$PRESSO.est,
                                     sei = subtypes_noTNBC$PRESSO.se)
  
  # Save results in I2_results
  I2_results[1,11] = round(meta_analysis_noTNBC.presso$I2, digits = 1)
  
  return(I2_results)
  
}

# Height
heterogeneity_height <- heterogeneity(risk_factor = "Height")
# BMI
heterogeneity_bmi <- heterogeneity(risk_factor = "Body mass index (kg/m2)")
# T2D
heterogeneity_t2d <- heterogeneity(risk_factor = "Type 2 diabetes")
# Age at menarche
heterogeneity_menarche <- heterogeneity(risk_factor = "Age at menarche (years)")
# Age at menopause
heterogeneity_menopause <- heterogeneity(risk_factor = "Age at menopause (years)")
# Percent density
heterogeneity_density <- heterogeneity(risk_factor = "Breast density (%)")
# Alcohol consumption
heterogeneity_alcohol <- heterogeneity(risk_factor = "Alcohol consumption (drinks/week)")
# Smoking
heterogeneity_smoking <- heterogeneity(risk_factor = "Smoking (ever smoked regularly; yes/no)")
# Physical activity
heterogeneity_activity <- heterogeneity(risk_factor = "Overall activity")


# Combine all estimates
all_I2 <- rbind(heterogeneity_height,
                heterogeneity_bmi,
                heterogeneity_t2d,
                heterogeneity_menarche,
                heterogeneity_menopause,
                heterogeneity_density,
                heterogeneity_alcohol,
                heterogeneity_smoking,
                heterogeneity_activity)

all_I2

# Remove individual I2 estimates
rm(heterogeneity_activity, heterogeneity_alcohol, heterogeneity_bmi, heterogeneity_density, heterogeneity_height,
   heterogeneity_menarche, heterogeneity_menopause, heterogeneity_smoking, heterogeneity_t2d)



## the end of the script ##
