### The aim of this Rscript is to perform the MR-PRESSO analyses for all nine risk factors and for three female-specific risk factors ###

# Set working directory
setwd("/DATA/users/r.verdiesen/MR_BCAC/")

# Load required libraries
library(dplyr)
library(parallel)
library(tictoc)

# Load summary level data for Mendelian randomization analyses (includes both exposure and outcome data)
load(file = "manuscript/allMRdata_finalversion.Rdata")


# Make list comprising all six outcome dataframes
outcomes <- list(overallbc_data, lumA_data, lumBHER2neg_data, lumB_data, HER2enr_data, TNBC_data)

names(outcomes) <- c("overallbc_data", "lumA_data", "lumBHER2neg_data", "lumB_data", "HER2enr_data", "TNBC_data")
# checked: names list are correct!


#################################### height ####################################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_height <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = height_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_height,
            sure = T)


# Create dataframe to store relevant estimates in
results_height_presso <- as.data.frame(matrix(nrow = 6,
                                              ncol = 6))

names(results_height_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_height_presso[i, 1] = "Height"
  
  results_height_presso[i, 2] = names(presso_height[i])
  
  results_height_presso[i, 3] = presso_height[[i]]$`Main MR results`[2,3]
  
  results_height_presso[i, 4] = presso_height[[i]]$`Main MR results`[2,4]
  
  results_height_presso[i, 5] = presso_height[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_height_presso[i, 6] = length(presso_height[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_height_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_height, results_height_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_height.RData")


# Quit R
quit(save = "no")


########################## height, female-specific data ########################

# Load female-specific height data into R
load("sumstat_riskfactors/female_betas/height_mrobject_women.RData")


# Clump female-specific data based on sex-combined data
height_mrobject_clumped <- TwoSampleMR::clump_data(dat = height_mrobject_def,
                                                   clump_kb = 1000,
                                                   clump_r2 = 0.001,
                                                   clump_p1 = 4.99e-08,
                                                   clump_p2 = 4.99e-08,
                                                   pop = "EUR")

height_mrobject_women_clumped <- height_mrobject_women %>%
  filter(SNP %in% height_mrobject_clumped$SNP)


# Remove in-between-objects from R
rm(height_mrobject_clumped,
   height_mrobject_def,
   height_mrobject_women)



# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_height_female <- mclapply(outcomes, function(outcome_data) {
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = height_mrobject_women_clumped,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_height_female,
            sure = T)


# Create dataframe to store relevant estimates in
results_height_female_presso <- as.data.frame(matrix(nrow = 6,
                                                     ncol = 6))

names(results_height_female_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_height_female_presso[i, 1] = "Height; female-specific data"
  
  results_height_female_presso[i, 2] = names(presso_height_female[i])
  
  results_height_female_presso[i, 3] = presso_height_female[[i]]$`Main MR results`[2,3]
  
  results_height_female_presso[i, 4] = presso_height_female[[i]]$`Main MR results`[2,4]
  
  results_height_female_presso[i, 5] = presso_height_female[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_height_female_presso[i, 6] = length(presso_height_female[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_height_female_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_height_female, results_height_female_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_height_femalespecific.RData")


# Quit R
quit(save = "no")


################################# BMI ##########################################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_bmi <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = bmi_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_bmi,
            sure = T)


# Create dataframe to store relevant estimates in
results_bmi_presso <- as.data.frame(matrix(nrow = 6,
                                           ncol = 6))

names(results_bmi_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_bmi_presso[i, 1] = "Body mass index (kg/m2)"
  
  results_bmi_presso[i, 2] = names(presso_bmi[i])
  
  results_bmi_presso[i, 3] = presso_bmi[[i]]$`Main MR results`[2,3]
  
  results_bmi_presso[i, 4] = presso_bmi[[i]]$`Main MR results`[2,4]
  
  results_bmi_presso[i, 5] = presso_bmi[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_bmi_presso[i, 6] = length(presso_bmi[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_bmi_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_bmi, results_bmi_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_bmi.RData")


# Quit R
quit(save = "no")


######################### BMI, female-specific data ############################

# Load female-specific bmi data into R
load("sumstat_riskfactors/female_betas/bmi_mrobject_women.RData")


# Clump female-specific data based on sex-combined data
bmi_mrobject_clumped <- TwoSampleMR::clump_data(dat = bmi_mrobject_def,
                                                clump_kb = 1000,
                                                clump_r2 = 0.001,
                                                clump_p1 = 4.99e-08,
                                                clump_p2 = 4.99e-08,
                                                pop = "EUR")

bmi_mrobject_women_clumped <- bmi_mrobject_women %>%
  filter(SNP %in% bmi_mrobject_clumped$SNP)


# Remove in-between-objects from R
rm(bmi_mrobject_clumped,
   bmi_mrobject_def,
   bmi_mrobject_women)



# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_bmi_female <- mclapply(outcomes, function(outcome_data) {
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = bmi_mrobject_women_clumped,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_bmi_female,
            sure = T)


# Create dataframe to store relevant estimates in
results_bmi_female_presso <- as.data.frame(matrix(nrow = 6,
                                                  ncol = 6))

names(results_bmi_female_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_bmi_female_presso[i, 1] = "BMI; female-specific data"
  
  results_bmi_female_presso[i, 2] = names(presso_bmi_female[i])
  
  results_bmi_female_presso[i, 3] = presso_bmi_female[[i]]$`Main MR results`[2,3]
  
  results_bmi_female_presso[i, 4] = presso_bmi_female[[i]]$`Main MR results`[2,4]
  
  results_bmi_female_presso[i, 5] = presso_bmi_female[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_bmi_female_presso[i, 6] = length(presso_bmi_female[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_bmi_female_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_bmi_female, results_bmi_female_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_bmi_femalespecific.RData")


# Quit R
quit(save = "no")


################################# T2D ##########################################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_t2d <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = t2d_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_t2d,
            sure = T)


# Create dataframe to store relevant estimates in
results_t2d_presso <- as.data.frame(matrix(nrow = 6,
                                           ncol = 6))

names(results_t2d_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_t2d_presso[i, 1] = "Type 2 diabetes"
  
  results_t2d_presso[i, 2] = names(presso_t2d[i])
  
  results_t2d_presso[i, 3] = presso_t2d[[i]]$`Main MR results`[2,3]
  
  results_t2d_presso[i, 4] = presso_t2d[[i]]$`Main MR results`[2,4]
  
  results_t2d_presso[i, 5] = presso_t2d[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_t2d_presso[i, 6] = length(presso_t2d[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_t2d_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_t2d, results_t2d_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_t2d.RData")


# Quit R
quit(save = "no")


############################# age at menarche ##################################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_menarche <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = menarche_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_menarche,
            sure = T)


# Create dataframe to store relevant estimates in
results_menarche_presso <- as.data.frame(matrix(nrow = 6,
                                                ncol = 6))

names(results_menarche_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_menarche_presso[i, 1] = "Age at menarche (years)"
  
  results_menarche_presso[i, 2] = names(presso_menarche[i])
  
  results_menarche_presso[i, 3] = presso_menarche[[i]]$`Main MR results`[2,3]
  
  results_menarche_presso[i, 4] = presso_menarche[[i]]$`Main MR results`[2,4]
  
  results_menarche_presso[i, 5] = presso_menarche[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_menarche_presso[i, 6] = length(presso_menarche[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_menarche_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_menarche, results_menarche_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_menarche.RData")


# Quit R
quit(save = "no")


########################### age at menopause ###################################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_menopause <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = menopause_mrobject_def[!is.na(menopause_mrobject_def$beta.exposure) ,], # otherwise the 2 NAs result in an error
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_menopause,
            sure = T)


# Create dataframe to store relevant estimates in
results_menopause_presso <- as.data.frame(matrix(nrow = 6,
                                                 ncol = 6))

names(results_menopause_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_menopause_presso[i, 1] = "Age at menopause (years)"
  
  results_menopause_presso[i, 2] = names(presso_menopause[i])
  
  results_menopause_presso[i, 3] = presso_menopause[[i]]$`Main MR results`[2,3]
  
  results_menopause_presso[i, 4] = presso_menopause[[i]]$`Main MR results`[2,4]
  
  results_menopause_presso[i, 5] = presso_menopause[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_menopause_presso[i, 6] = length(presso_menopause[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_menopause_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_menopause, results_menopause_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_menopause.RData")


# Quit R
quit(save = "no")


########################### percent breast density #############################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_density <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = density_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_density,
            sure = T)


# Create dataframe to store relevant estimates in
results_density_presso <- as.data.frame(matrix(nrow = 6,
                                               ncol = 6))

names(results_density_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_density_presso[i, 1] = "Breast density (%)"
  
  results_density_presso[i, 2] = names(presso_density[i])
  
  results_density_presso[i, 3] = presso_density[[i]]$`Main MR results`[2,3]
  
  results_density_presso[i, 4] = presso_density[[i]]$`Main MR results`[2,4]
  
  results_density_presso[i, 5] = presso_density[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_density_presso[i, 6] = length(presso_density[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_density_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_density, results_density_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_density.RData")


# Quit R
quit(save = "no")


################################ alcohol consumption ###########################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_alcohol <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = alcohol_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_alcohol,
            sure = T)


# Create dataframe to store relevant estimates in
results_alcohol_presso <- as.data.frame(matrix(nrow = 6,
                                               ncol = 6))

names(results_alcohol_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_alcohol_presso[i, 1] = "Alcohol consumption (drinks/week)"
  
  results_alcohol_presso[i, 2] = names(presso_alcohol[i])
  
  results_alcohol_presso[i, 3] = presso_alcohol[[i]]$`Main MR results`[2,3]
  
  results_alcohol_presso[i, 4] = presso_alcohol[[i]]$`Main MR results`[2,4]
  
  results_alcohol_presso[i, 5] = presso_alcohol[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_alcohol_presso[i, 6] = length(presso_alcohol[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_alcohol_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_alcohol, results_alcohol_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_alcohol.RData")


# Quit R
quit(save = "no")


################################ smoking #######################################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_smoking <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = smoking_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_smoking,
            sure = T)


# Create dataframe to store relevant estimates in
results_smoking_presso <- as.data.frame(matrix(nrow = 6,
                                               ncol = 6))

names(results_smoking_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_smoking_presso[i, 1] = "Smoking (ever smoked regularly; yes/no)"
  
  results_smoking_presso[i, 2] = names(presso_smoking[i])
  
  results_smoking_presso[i, 3] = presso_smoking[[i]]$`Main MR results`[2,3]
  
  results_smoking_presso[i, 4] = presso_smoking[[i]]$`Main MR results`[2,4]
  
  results_smoking_presso[i, 5] = presso_smoking[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_smoking_presso[i, 6] = length(presso_smoking[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_smoking_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_smoking, results_smoking_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_smoking.RData")


# Quit R
quit(save = "no")


################################# phisical activity ############################

# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_activity <- mclapply(outcomes, function(outcome_data) {
  
  # Clump genetic instruments for risk factor
  clumped_data <- TwoSampleMR::clump_data(dat = activity_mrobject_def,
                                          clump_kb = 1000,
                                          clump_r2 = 0.001,
                                          clump_p1 = 4.99e-08,
                                          clump_p2 = 4.99e-08,
                                          pop = "EUR")
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = clumped_data,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_activity,
            sure = T)


# Create dataframe to store relevant estimates in
results_activity_presso <- as.data.frame(matrix(nrow = 6,
                                                ncol = 6))

names(results_activity_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_activity_presso[i, 1] = "Overall activity"
  
  results_activity_presso[i, 2] = names(presso_activity[i])
  
  results_activity_presso[i, 3] = presso_activity[[i]]$`Main MR results`[2,3]
  
  results_activity_presso[i, 4] = presso_activity[[i]]$`Main MR results`[2,4]
  
  results_activity_presso[i, 5] = presso_activity[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_activity_presso[i, 6] = length(presso_activity[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_activity_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_activity, results_activity_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_activity.RData")


# Quit R
quit(save = "no")


##################### physical activity, female-specific data ##################

# Load female-specific activity data into R
load("sumstat_riskfactors/female_betas/activity_mrobject_women.RData")


# Clump female-specific data based on sex-combined data
activity_mrobject_clumped <- TwoSampleMR::clump_data(dat = activity_mrobject_def,
                                                     clump_kb = 1000,
                                                     clump_r2 = 0.001,
                                                     clump_p1 = 4.99e-08,
                                                     clump_p2 = 4.99e-08,
                                                     pop = "EUR")

activity_mrobject_women_clumped <- activity_mrobject_women %>%
  filter(SNP %in% activity_mrobject_clumped$SNP)


# Remove in-between-objects from R
rm(activity_mrobject_clumped,
   activity_mrobject_def,
   activity_mrobject_women)



# Run MR PRESSO pipeline function for each outcome in parallel using mclapply function

# Start time benchmark
tic()

# Make sure that estimates will be reproducible
RNGkind("L'Ecuyer-CMRG")

# mclapply function
presso_activity_female <- mclapply(outcomes, function(outcome_data) {
  
  
  # Harmonize risk factor-outcome data
  harmon_data <- TwoSampleMR::harmonise_data(exposure_dat = activity_mrobject_women_clumped,
                                             outcome_dat = outcome_data,
                                             action = 2)
  
  
  # Exclude palindromic SNPs with intermediate allele frequencies before running MR-PRESSO
  harmon_data_2 <- harmon_data %>%
    filter(mr_keep == "TRUE")
  
  
  # Perform MR-PRESSO
  MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                      BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome",
                      SdExposure = "se.exposure",
                      data = harmon_data_2,
                      OUTLIERtest = T,
                      DISTORTIONtest = T,
                      SignifThreshold = 0.05,
                      NbDistribution = 60000,
                      seed = 43890)
  
}, mc.cores = 6)


# End time benchmark
toc()


# Clean up R environment before saving relevant objects
gdata::keep(presso_activity_female,
            sure = T)


# Create dataframe to store relevant estimates in
results_activity_female_presso <- as.data.frame(matrix(nrow = 6,
                                                       ncol = 6))

names(results_activity_female_presso) <- c("Risk_factor", "Outcome", "PRESSO.est", "PRESSO.se", "PRESSO.pleio.p", "PRESSO.n.outliers")


# Use for loop to store relevant estimates into results object
for (i in 1:6){
  
  results_activity_female_presso[i, 1] = "Overall physical activity; female-specific data"
  
  results_activity_female_presso[i, 2] = names(presso_activity_female[i])
  
  results_activity_female_presso[i, 3] = presso_activity_female[[i]]$`Main MR results`[2,3]
  
  results_activity_female_presso[i, 4] = presso_activity_female[[i]]$`Main MR results`[2,4]
  
  results_activity_female_presso[i, 5] = presso_activity_female[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  
  results_activity_female_presso[i, 6] = length(presso_activity_female[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  
}

results_activity_female_presso # checked; correct


# Clean up R environment before saving MR-PRESSO objects and results
gdata::keep(presso_activity_female, results_activity_female_presso,
            sure = T)


# Save R environment
save.image(file = "manuscript/Robjects_mrpresso_activity_femalespecific.RData")


# Quit R
quit(save = "no")


## the end of the script ##
