## The aim of this Rscript is to create tables in the main manuscript, and in the supplemental file ##

############################# manuscript table 2 ###############################
# Subtype-specific causal effect estimates per unit increase for all nine breast cancer risk factors across primary and secondary MR analyses.
# ORs correspond to a 1 standard deviation increase for all risk factors, except for age at menarche and age at menopause. ORs for these two risk factors correspond to a 1 year increase.

# Bind all elements for table 2 into 1 dataframe
all_riskfactors_table2 <- rbind(height_all_table2,
                                bmi_all_table2,
                                t2d_all_table2,
                                menarche_all_table2,
                                menopause_all_table2,
                                density_all_table2,
                                alcohol_all_table2,
                                smoking_all_table2,
                                activity_all_table2)

# note: the inpute dara frames used for making table2 are the rusults of primary and secondary analysis for each risk factor. 
# To know how they were created, check the Rscript "MR_Primary_And_Secondary_Analysis".


########################## supplemental table 3 ################################
# Calculated F-statistics for each included risk factor.

# Make function for calculation F statistics
# The formula for this is: F statistic = R^2/(1-R^2)*sample size/number of variants
# The r2 and nsnps values were extracted from the original GWAS papers
# fstat_min is the F-stat corresponding to analyses with HER-enrinched tumors as outcome
# fstat_max is the F-stat corresponding to analyses with Luminal A-like tumors as outcome


fstat <- function(r2, nsnps){
  
  fstat_min <- r2/(1-r2)*94361/nsnps
  
  fstat_max <- r2/(1-r2)*136730/nsnps
  
  return(paste(round(fstat_min, digits = 1), round(fstat_max, digits = 1), sep = ", "))
  
}


# Calculate risk factor specific F-stats

# Height
fstat_height <- fstat(r2 = 0.246, nsnps = 3290)

# BMI
fstat_bmi <- fstat(r2 = 0.06, nsnps = 941)

# T2D
# Cannot be calculated because r2 was not reported in GWAS

# Age at menarche
fstat_menarche <- fstat(r2 = 0.074, nsnps = 389)

# Age at menopause
fstat_menopause <- fstat(r2 = 0.101, nsnps = 290)

# Percent breast density
fstat_density <- fstat(r2 = 0.087, nsnps = 20)

# Alcohol consumption
fstat_alcohol <- fstat(r2 = 0.0019, nsnps = 99)

# Smoking
fstat_smoking <- fstat(r2 = 0.023, nsnps = 378)

# Activity
fstat_activity <- fstat(r2 = 0.0006, nsnps = 3)


# The results of these calculations are presented in Supplemental Table 3

rm(list = ls())



############################### supplemental table 4 ###########################
# Results of primary and secondary MR analyses with female-specific estimates for height, BMI and physical activity

### Height ###

# Load PCA results for female-specific analyses
load("manuscript/data_height_femalespecific.Rdata")

height_pca_all_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, IVW.est, CI.low, CI.high) %>%
  mutate(OR = round(exp(IVW.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(CI.low), digits = 2), round(exp(CI.high), digits = 2), sep = ", "))


# Load results from IVW analyses with clumped variants into R
load("manuscript/data_height_femalespecific_clumped_dd13122022.Rdata")

# IVW with clumped variants
height_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, IVW.est, CI.low, CI.high) %>%
  mutate(OR = round(exp(IVW.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(CI.low), digits = 2), round(exp(CI.high), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

# MR-Egger
height_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Egger.est, Egger.se, P.pleiotr) %>%
  mutate(OR = round(exp(Egger.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Egger.est - (1.96*Egger.se)), digits = 2), round(exp(Egger.est + (1.96*Egger.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR, P.pleiotr)

# Weighted median
height_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Median.est, Median.se) %>%
  mutate(OR = round(exp(Median.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Median.est - (1.96*Median.se)), digits = 2), round(exp(Median.est + (1.96*Median.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

rm(list = ls())

### BMI ###

# Load PCA results for female-specific analyses
load("manuscript/data_bmi_femalespecific.Rdata")

bmi_pca_all_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, IVW.est, CI.low, CI.high) %>%
  mutate(OR = round(exp(IVW.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(CI.low), digits = 2), round(exp(CI.high), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

rm(bmi_pca_all_female)

# Load results from IVW analyses with clumped variants into R
load("manuscript/data_bmi_femalespecific_clumped_dd13122022.Rdata")

# IVW
bmi_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, IVW.est, CI.low, CI.high) %>%
  mutate(OR = round(exp(IVW.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(CI.low), digits = 2), round(exp(CI.high), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

# MR-Egger
bmi_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Egger.est, Egger.se, P.pleiotr) %>%
  mutate(OR = round(exp(Egger.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Egger.est - (1.96*Egger.se)), digits = 2), round(exp(Egger.est + (1.96*Egger.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR, P.pleiotr)

# Weighted median
bmi_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Median.est, Median.se) %>%
  mutate(OR = round(exp(Median.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Median.est - (1.96*Median.se)), digits = 2), round(exp(Median.est + (1.96*Median.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

# Weighted mode
bmi_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Mode.est, Mode.se) %>%
  mutate(OR = round(exp(Mode.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Mode.est - (1.96*Mode.se)), digits = 2), round(exp(Mode.est + (1.96*Mode.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

load("manuscript/Robjects_mrpresso_bmi_femalespecific.RData")

results_bmi_female_presso %>%
  filter(Outcome != "overallbc_data") %>%
  select(Risk_factor, Outcome, PRESSO.est, PRESSO.se) %>%
  select(Risk_factor, Outcome, OR, PRESSO.CI)

## 95%CI luminal is incorrect!!

rm(list = ls())

### Physical activity ###

# Load PCA results for female-specific analyses
load("manuscript/data_activity_femalespecific.Rdata")

activity_pca_all_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, IVW.est, CI.low, CI.high) %>%
  mutate(OR = round(exp(IVW.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(CI.low), digits = 2), round(exp(CI.high), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

rm(activity_pca_all_female)


# Load results from IVW analyses with clumped variants into R
load("manuscript/data_activity_femalespecific_clumped_dd13122022.Rdata")

# IVW
activity_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, IVW.est, CI.low, CI.high) %>%
  mutate(OR = round(exp(IVW.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(CI.low), digits = 2), round(exp(CI.high), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

# MR-Egger
activity_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Egger.est, Egger.se, P.pleiotr) %>%
  mutate(OR = round(exp(Egger.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Egger.est - (1.96*Egger.se)), digits = 2), round(exp(Egger.est + (1.96*Egger.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR, P.pleiotr)

# Weighted median
activity_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Median.est, Median.se) %>%
  mutate(OR = round(exp(Median.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Median.est - (1.96*Median.se)), digits = 2), round(exp(Median.est + (1.96*Median.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

# Weighted mode
activity_all_clumped_female %>%
  filter(Outcome != "All BCAC breast cancer cases") %>%
  select(`Risk factor`, Outcome, Mode.est, Mode.se) %>%
  mutate(OR = round(exp(Mode.est), digits = 2)) %>%
  mutate(CI.OR = paste(round(exp(Mode.est - (1.96*Mode.se)), digits = 2), round(exp(Mode.est + (1.96*Mode.se)), digits = 2), sep = ", ")) %>%
  select(`Risk factor`, Outcome, OR, CI.OR)

rm(list = ls())

# MR PRESSO
load("manuscript/Robjects_mrpresso_activity_femalespecific.RData")

results_activity_female_presso %>%
  filter(Outcome != "overallbc_data") %>%
  select(Risk_factor, Outcome, PRESSO.est, PRESSO.se) %>%
  mutate(OR = round(exp(PRESSO.est), digits = 2)) %>%
  mutate(PRESSO.CI = paste(round(exp(PRESSO.est - (1.96*PRESSO.se)), digits = 2), round(exp(PRESSO.est + (1.96*PRESSO.se)), digits = 2), sep = ", ")) %>%
  select(Risk_factor, Outcome, OR, PRESSO.CI)

## 95%CI luminal is incorrect!!
exp(-0.1116748 - (1.96*0.02274853)); exp(-0.1116748 + (1.96*0.02274853))


## the end of the script ##
