# The aim of this Rscript is to calculate % of heterogeneity across subtype-specific estimates

load("/DATA/projects/mr_bcac/manuscript/data_fig1_dd02022022.Rdata")

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
