## The aim of this Rscript is to create figures in the main manuscript, and in the supplemental file ##

############################# manuscript figure 1 ##############################
# Causal breast cancer subtype-specific effect estimates per unit increase for nine established breast cancer risk factors. 

# Load required libraries 
library(metafor)
library(ggpubr)
library(dplyr)

# Load data with estimates for figure 1 into R
load("/DATA/projects/mr_bcac/manuscript/data_fig1_dd02022022.Rdata")

# Combine data for all risk factors into one dataframe
all_riskfactors_fig1 = rbind(height_all_fig1,
                             bmi_all_fig1,
                             t2d_all_fig1,
                             menarche_all_fig1,
                             menopause_all_fig1,
                             density_all_fig1,
                             alcohol_all_fig1,
                             smoking_all_fig1,
                             activity_all_fig1)


# Edit the name of overall activity to physical activity
activity_all_fig1 <- activity_all_fig1 %>%
  mutate(`Risk factor` = replace(`Risk factor`, `Risk factor` == "Overall activity", "Physical activity"))


# Subset data to results including PCs only
fig1_pcs = subset(all_riskfactors_fig1,
                  all_riskfactors_fig1$Analysis == "IVW including PCs")


# Load data with IVW estimates clumped variants into R
load("/DATA/projects/mr_bcac/manuscript/data_table2_dd27062022.Rdata")

# Edit the name of overall activity to physical activity
all_riskfactors_table2 <- all_riskfactors_table2 %>%
  mutate(`Risk factor` = replace(`Risk factor`, `Risk factor` == "Overall activity", "Physical activity"))


head(all_riskfactors_table2)

# Subset all_riskfactors_table2 to IVW results only
ivw_clumped <- all_riskfactors_table2[, c(1:7)]

ivw_clumped$se.est <- NA
ivw_clumped$Analysis <- "IVW clumped IVs"
ivw_clumped$order <- seq(1:54)
ivw_clumped$N <- NA


# Make an additional variable containing the number of cases for each subtype
fig1_pcs$N <- ifelse(fig1_pcs$Outcome == "All BCAC breast cancer cases", "133,384",
                     ifelse(fig1_pcs$Outcome == "Luminal A-like", "45,253",
                            ifelse(fig1_pcs$Outcome == "Luminal B/HER2-negative like", "6,350",
                                   ifelse(fig1_pcs$Outcome == "Luminal B-like", "6,427",
                                          ifelse(fig1_pcs$Outcome == "HER2-enriched-like", "2,884", "8,602")))))

# Make an additional variable to order the rows in the forest plot on
fig1_pcs$order <- seq(1:54)

# Edit the name of overall activity to physical activity
fig1_pcs <- fig1_pcs %>%
  mutate(`Risk factor` = replace(`Risk factor`, `Risk factor` == "Overall activity", "Physical activity")) %>%
  rename(IVW.se = se.est, se.est = p.hetr.snps)


# Merge data to fig1_pcs dataframe
all_ivw_est <- rbind(fig1_pcs,
                     ivw_clumped)


all_ivw_est_fig1 <- all_ivw_est[order(all_ivw_est$order), ]


# Make order_plot variable
all_ivw_est_fig1$order_plot <- seq(1:108)

all_ivw_est_fig1$order <- NULL

fig1_pcs$order <- rep(c(1:6), 9)

anthropometric <- fig1_pcs %>%
  filter(`Risk factor` == "Height" |  
           `Risk factor` == "Body mass index (kg/m2)" | 
           `Risk factor` == "Type 2 diabetes")

anthropometric$rf <- factor(anthropometric$`Risk factor`, levels = c("Height", "Body mass index (kg/m2)", "Type 2 diabetes"))

forest1 <- ggplot(anthropometric, aes(y = reorder(Outcome, -order), x = exp(IVW.est))) +
  geom_point(size = 2) +
  coord_trans(x = "log10") +
  #scale_colour_manual("Molecular subtype", values = c(inauguration("inauguration_2021_bernie")[1:6])) +
  #, breaks = c("TN", "LumB", "LumA", "HER2")) +
  facet_wrap(~ rf, ncol = 3, scales = "fixed", labeller = labeller(rf = 
                                                                     c("Height" = "Height", 
                                                                       "Body mass index (kg/m2)" = "Body mass index (kg/m2)", 
                                                                       "Type 2 diabetes" = "Type 2 diabetes"))) +
  theme_bw() +
  theme(legend.position="none") +
  geom_segment(aes(x = exp(CI.low), xend = exp(CI.high), yend = Outcome), size = 1) +
  geom_vline(lty = 2, aes(xintercept = 1), colour = "grey") +
  xlab("Causal odds ratio") +
  ylab("") +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold"))

forest1


reproductive <- fig1_pcs %>%
  filter(`Risk factor` == "Age at menarche (years)" |
           `Risk factor` == "Age at menopause (years)" |
           `Risk factor` == "Breast density (%)")

forest2 <- ggplot(reproductive, aes(y = reorder(Outcome, -order), x = exp(IVW.est))) +
  geom_point(size = 2) +
  coord_trans(x = "log10") +
  #scale_colour_manual("Molecular subtype", values = c(inauguration("inauguration_2021_bernie")[1:6])) +
  #, breaks = c("TN", "LumB", "LumA", "HER2")) +
  facet_wrap(~ `Risk factor`, ncol = 3, scales = "fixed", labeller = labeller(`Risk factor` = 
                                                                                c("Age at menarche (years)" = "Age at menarche (years)", 
                                                                                  "Age at menopause (years)" = "Age at menopause (years)", 
                                                                                  "Breast density (%)" = "Breast density (%)"))) +
  theme_bw() +
  theme(legend.position="none") +
  geom_segment(aes(x = exp(CI.low), xend = exp(CI.high), yend = Outcome), size = 1) +
  geom_vline(lty = 2, aes(xintercept = 1), colour = "grey") +
  xlab("Causal odds ratio") +
  ylab("") +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold"))

forest2


lifestyle <- fig1_pcs %>%
  filter(`Risk factor` == "Alcohol consumption (drinks/week)" |
           `Risk factor` == "Smoking (ever smoked regularly; yes/no)" |
           `Risk factor` == "Physical activity")

lifestyle$rf <- factor(lifestyle$`Risk factor`, levels = c("Alcohol consumption (drinks/week)", "Smoking (ever smoked regularly; yes/no)", "Physical activity"))

forest3 <- ggplot(lifestyle, aes(y = reorder(Outcome, -order), x = exp(IVW.est))) +
  geom_point(size = 2) +
  coord_trans(x = "log10") +
  #scale_colour_manual("Molecular subtype", values = c(inauguration("inauguration_2021_bernie")[1:6])) +
  #, breaks = c("TN", "LumB", "LumA", "HER2")) +
  facet_wrap(~ rf, ncol = 3, scales = "fixed", labeller = labeller(rf = 
                                                                     c("Alcohol consumption (drinks/week)" = "Alcohol consumption (drinks/week)", 
                                                                       "Smoking (ever smoked regularly; yes/no)" = "Smoking (ever smoked regularly; yes/no)", 
                                                                       "Physical activity" = "Physical activity"))) +
  theme_bw() +
  theme(legend.position="none") +
  geom_segment(aes(x = exp(CI.low), xend = exp(CI.high), yend = Outcome), size = 1) +
  geom_vline(lty = 2, aes(xintercept = 1), colour = "grey") +
  xlab("Causal odds ratio") +
  ylab("") +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold"))

forest3


# merge three plots into 1 using ggarange = figure 1
fig1 <- ggarrange(plotlist = list(forest1, forest2, forest3),
                  nrow = 3, ncol = 1)



############################# manuscript figure 2 ##############################
# Overview of evidence for (subtype-specific) causal effects per increasing unit of the risk factor. 
# Define categories in each depth: 

# Depth 1: Strength of causal evidence based on combination of primary and secondary MR analysis. (shown by size)
#           1. robust: if all MR estimates reached p < 0.05 and all estimates were concordant in direction. 
#           2. probable: if at least one method reached p < 0.05 and all effect estimates were concordant in direction. 
#           3. suggestive: if at least one method reached p < 0.05 but effect estimates were not concordant in direction. 
#           4. insufficient: Empty cells, evidence for a causal association was insufficient. 

# Depth 2: Direction of causal effect (shown by color)
#           1. -1: risk decreasing
#           2.  0: inconclusive 
#           3. +1: risk increasing

# Depth 3: betas extracted from summary statistics GWAS for MR analysis
#           1: female-specific betas (solid color)
#           2. sex-combined betas (fade color)

# Create an array includes all nine risk factors, all are double named mr_overview.rds

# Define dimensions
rows <- c("Height", "Height", 
          "BMI", "BMI", 
          "T2D", "T2D", 
          "Age at menarche", "Age at menarche", 
          "Age at menopause","Age at menopause", 
          "Percent breast density", "Percent breast density",
          "Alcohol consumption", "Alcohol consumption", 
          "Smoking behaviour", "Smoking behaviour",
          "Physical activity", "Physical activity")

cols <- c("Luminal A", "Luminal B/HER2-negative", "Luminal B", "HER2-enriched", "Triple negative")

# Create the data matrix with dimension names
mr_overview <- array(18, dim = c(length(rows), length(cols), 3))  # 3 depths

# Assign values to the matrix
# row 1 includes Height, female specific betas
mr_overview[1, 1, ] <- c(2, +1, 1) # c() represent three depths in the array
mr_overview[1, 2, ] <- c(2, +1, 1)
mr_overview[1, 3, ] <- c(1, +1, 1)
mr_overview[1, 4, ] <- c(2, +1, 1)
mr_overview[1, 5, ] <- c(2, +1, 1)

# row 2 includes Height, sex-combined betas
mr_overview[2, 1, ] <- c(1, +1, 2) 
mr_overview[2, 2, ] <- c(2, +1, 2)
mr_overview[2, 3, ] <- c(2, +1, 2)
mr_overview[2, 4, ] <- c(4, 9, 2)
mr_overview[2, 5, ] <- c(4, 9, 2)

# row 3 includes BMI, female-specific betas
mr_overview[3, 1, ] <- c(2, -1, 1) 
mr_overview[3, 2, ] <- c(2, -1, 1)
mr_overview[3, 3, ] <- c(2, -1, 1)
mr_overview[3, 4, ] <- c(2, -1, 1)
mr_overview[3, 5, ] <- c(2, -1, 1)

# row 4 includes BMI, sex-combined betas
mr_overview[4, 1, ] <- c(4, 9, 2) 
mr_overview[4, 2, ] <- c(4, 9, 2)
mr_overview[4, 3, ] <- c(3, 0, 2)
mr_overview[4, 4, ] <- c(3, 0, 2)
mr_overview[4, 5, ] <- c(4, 9, 2)

# row 5 includes T2D female-specific betas
mr_overview[5, 1, ] <- c(-9, -9, -9) 
mr_overview[5, 2, ] <- c(-9, -9, -9)
mr_overview[5, 3, ] <- c(-9, -9, -9)
mr_overview[5, 4, ] <- c(-9, -9, -9)
mr_overview[5, 5, ] <- c(-9, -9, -9)

# row 6 includes T2D sex-combined betas
mr_overview[6, 1, ] <- c(3, 0, 2) 
mr_overview[6, 2, ] <- c(4, 9, 2)
mr_overview[6, 3, ] <- c(4, 9, 2)
mr_overview[6, 4, ] <- c(2, -1, 2)
mr_overview[6, 5, ] <- c(3, 0, 2)

# row 7 includes Age at menarche, female-specific betas 
mr_overview[7, 1, ] <- c(2, -1, 1) 
mr_overview[7, 2, ] <- c(3, 0, 1)
mr_overview[7, 3, ] <- c(4, 9, 1)
mr_overview[7, 4, ] <- c(4, 9, 1)
mr_overview[7, 5, ] <- c(3, 0, 1)

# row 8 includes Age at menarche, sex-combined betas 
mr_overview[8, 1, ] <- c(-9, -9, -9) 
mr_overview[8, 2, ] <- c(-9, -9, -9)
mr_overview[8, 3, ] <- c(-9, -9, -9)
mr_overview[8, 4, ] <- c(-9, -9, -9)
mr_overview[8, 5, ] <- c(-9, -9, -9)

# row 9 includes Age at menopause, female-specific betas
mr_overview[9, 1, ] <- c(1, +1, 1) 
mr_overview[9, 2, ] <- c(2, +1, 1)
mr_overview[9, 3, ] <- c(2, +1, 1)
mr_overview[9, 4, ] <- c(2, +1, 1)
mr_overview[9, 5, ] <- c(4, 9, 1)

# row 10 includes Age at menopause, sex-combined betas
mr_overview[10, 1, ] <- c(-9, -9, -9) 
mr_overview[10, 2, ] <- c(-9, -9, -9)
mr_overview[10, 3, ] <- c(-9, -9, -9)
mr_overview[10, 4, ] <- c(-9, -9, -9)
mr_overview[10, 5, ] <- c(-9, -9, -9)

# row 11 includes Percent breast density, female-specific betas
mr_overview[11, 1, ] <- c(3, 0, 1) 
mr_overview[11, 2, ] <- c(3, 0, 1)
mr_overview[11, 3, ] <- c(4, 9, 1)
mr_overview[11, 4, ] <- c(4, 9, 1)
mr_overview[11, 5, ] <- c(4, 9, 1)

# row 12 includes Percent breast density, sex-combined betas
mr_overview[12, 1, ] <- c(-9, -9, -9) 
mr_overview[12, 2, ] <- c(-9, -9, -9)
mr_overview[12, 3, ] <- c(-9, -9, -9)
mr_overview[12, 4, ] <- c(-9, -9, -9)
mr_overview[12, 5, ] <- c(-9, -9, -9)

# row 13 includes Alcohol consumption, female-specific betas
mr_overview[13, 1, ] <- c(-9, -9, -9) 
mr_overview[13, 2, ] <- c(-9, -9, -9)
mr_overview[13, 3, ] <- c(-9, -9, -9)
mr_overview[13, 4, ] <- c(-9, -9, -9)
mr_overview[13, 5, ] <- c(-9, -9, -9)

# row 14 includes Alcohol consumption, sex-combined betas
mr_overview[14, 1, ] <- c(4, 9, 2) 
mr_overview[14, 2, ] <- c(4, 9, 2)
mr_overview[14, 3, ] <- c(2, -1, 2)
mr_overview[14, 4, ] <- c(4, 9, 2)
mr_overview[14, 5, ] <- c(4, 9, 2)

# row 15 includes Smoking behavior, female-specific betas
mr_overview[15, 1, ] <- c(-9, -9, -9) 
mr_overview[15, 2, ] <- c(-9, -9, -9)
mr_overview[15, 3, ] <- c(-9, -9, -9)
mr_overview[15, 4, ] <- c(-9, -9, -9)
mr_overview[15, 5, ] <- c(-9, -9, -9)

# row 16 includes Smoking behavior, sex-combined betas
mr_overview[16, 1, ] <- c(4, 9, 2)
mr_overview[16, 2, ] <- c(4, 9, 2)
mr_overview[16, 3, ] <- c(4, 9, 2)
mr_overview[16, 4, ] <- c(4, 9, 2)
mr_overview[16, 5, ] <- c(4, 9, 2)

# row 17 includes Physical activity, female-specific betas 
mr_overview[17, 1, ] <- c(2, -1, 1) 
mr_overview[17, 2, ] <- c(2, -1, 1)
mr_overview[17, 3, ] <- c(3, 0, 1)
mr_overview[17, 4, ] <- c(2, -1, 1)
mr_overview[17, 5, ] <- c(4, 9, 1)

# row 18 includes Physical activity, sex-combined betas
mr_overview[18, 1, ] <- c(2, -1, 2)
mr_overview[18, 2, ] <- c(3, 0, 2)
mr_overview[18, 3, ] <- c(4, 9, 2)
mr_overview[18, 4, ] <- c(4, 9, 2)
mr_overview[18, 5, ] <- c(4, 9, 2)

# Assign row and column names
dimnames(mr_overview) <- list(rows, cols, NULL)
print(mr_overview) # checked
saveRDS(mr_overview, "/DATA/users/m.shokouhi/projects/MR_BCsubtypes/mr_overview.rds")

# grid data
## colors
ColorNames  <- c('green', 'blue', 'red','black')

## griding space
grid_data <- array(9, dim = c(length(rows), length(cols), 5))  # 3 depths; 1:X, 2:Y, 3:color, 4:shape, 5:size

## X axis
x           <- seq(1, 2 * length(cols) - 1, by=2)
XTicks      <- x
XTickLabels <- c("Luminal A"                 , 
                 "Luminal B/ \nHER2-negative", 
                 "Luminal B"                 , 
                 "HER2-enriched"             , 
                 "Triple negative"           )

## Y axis
y           <- seq(length(rows)    , 1  , by=-1)
YTicks      <- seq(length(rows)-0.5, 1.5, by=-2)
YTickLabels <- c("Height"                , 
                 "BMI"                   , 
                 "T2D"                   , 
                 "Age at menarche"       , 
                 "Age at menopause"      , 
                 "Percent breast density", 
                 "Alcohol consumption"   , 
                 "Smoking behaviour"     , 
                 "Physical activity"     )


for (i in 1:length(rows)) {
  for (j in 1:length(cols)) {
    # X
    grid_data[i, j, 1] <- x[j]
    
    # Y
    grid_data[i, j, 2] <- y_1[i]
    
    # color
    if (as.numeric(mr_overview[i, j, 2])        == -1) {  # decreasing
      grid_data[i, j, 3] <- ColorNames[1]
    } else if (as.numeric(mr_overview[i, j, 2]) ==  0) {  # inconclusive
      grid_data[i, j, 3] <- ColorNames[2]
    } else if (as.numeric(mr_overview[i, j, 2]) ==  1) {  # increasing
      grid_data[i, j, 3] <- ColorNames[3]
    } else if (as.numeric(mr_overview[i, j, 2]) ==  9) {  # non-significant
      grid_data[i, j, 3] <- ColorNames[4]
    } else if (as.numeric(mr_overview[i, j, 2]) == -9) {  # no data available
      grid_data[i, j, 3] <- ColorNames[4]
    }
    
    # shape
    if (mr_overview[i, j, 3] == 1) { 
      grid_data[i, j, 4] <- 19 # filled "O"
    } else if (mr_overview[i, j, 3] == 2) { 
      grid_data[i, j, 4] <- 19 # filled "O"
    } else {  
      grid_data[i, j, 4] <- 4 # "X"
    }
    # size
    if (as.numeric(mr_overview[i, j, 1]) != -9) {
      grid_data[i, j, 5] <- 4 - as.numeric(mr_overview[i, j, 1])
    } else {
      grid_data[i, j, 5] <- 1
    }
  }
}


# Show female-specific results in the first row and sex-combined in the second row in Height, BMI, and activity 
par(mar = c(3, 7, 0.001, .5) + 3)           # for the paper size and the positions of the plot on the paper.
plot(as.numeric(grid_data[,,1])           , # x values
     as.numeric(grid_data[,,2])           , # y values
     xlab = ""                            ,
     ylab = ""                            ,
     col = grid_data[,,3]                 , # color
     pch =  as.numeric(grid_data[,,4])    , # shape
     cex  = as.numeric(grid_data[,,5])*1.5, # size
     axes = FALSE                         , # dont draw the axis
     xlim = c(0.0,9.5)                    , # limits of the X axis 
     ylim = c(1,length(rows))             ) # limits of the Y axis 

# horizontal solid lines
for (i in seq(length(rows)-1.5,1.5,by=-2)){
  lines(c(-0.5,10), c(i,i)  , col = "black", type = "l", lwd = 2.0, lty = 1)
}
# horizontal dashed lines
for (i in seq(length(rows)-0.5,1.5,by=-2)){
  lines(c(-0.5,10), c(i,i)  , col = "black", type = "l", lwd = 0.5, lty = 2)
}

# vertical lines
for (i in seq(8,0,by = -2)) {
  lines(c(i,i), c(0,length(rows) + 1), col = "black", type = "l", lwd = 2.0, lty = 1)
}

box()

axis(1, at = XTicks, labels = XTickLabels, las = 1, tcl = -0.02) # for x ticks
axis(2, at = YTicks, labels = YTickLabels, las = 1, tcl = -0.02) # for y ticks



############################### supplemental figure 2 ##########################
# Forrest plots with causal effects of nine risk factors on each hormone receptor breast cancer subtype across primary and secondary MR methods.

library(dplyr)
library(ggplot2)
# Load estimates robust methods into R
load("/DATA/projects/mr_bcac/manuscript/data_table2_dd27062022.Rdata") # includes all robust analysis results (except MR-PRESSO) for all risk factors (except female-specific)
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_height.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_bmi.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_t2d.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_menarche.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_menopause.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_density.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_smoking.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_alcohol.RData")
load("/DATA/projects/mr_bcac/manuscript/Robjects_mrpresso_activity.RData")
# Also add estimates for primary analysis
load("/DATA/projects/mr_bcac/manuscript/data_fig1_dd02022022.Rdata")


# Combine MR PRESSO results into 1 dataframe and merge to table 2 estimates
all_presso <- rbind(results_height_presso,
                    results_bmi_presso,
                    results_t2d_presso,
                    results_menarche_presso,
                    results_menopause_presso,
                    results_density_presso,
                    results_smoking_presso,
                    results_alcohol_presso,
                    results_activity_presso)

head(all_presso)


# Edit the name of overall activity to physical activity
activity_all_fig1 <- activity_all_fig1 %>%
  mutate(`Risk factor` = replace(`Risk factor`, `Risk factor` == "Overall activity", "Physical activity"))

all_presso <- all_presso %>%
  mutate(`Risk_factor` = replace(`Risk_factor`, `Risk_factor` == "Overall activity", "Physical activity"))

all_riskfactors_table2 <- all_riskfactors_table2 %>%
  mutate(`Risk factor` = replace(`Risk factor`, `Risk factor` == "Overall activity", "Physical activity"))


# Create outcome variable for merging of MR-PRESSO results to table 2 object
all_presso$Outcome_def <- ifelse(all_presso$Outcome == "overallbc_data", "All BCAC breast cancer cases",
                                 ifelse(all_presso$Outcome == "lumA_data", "Luminal A-like",
                                        ifelse(all_presso$Outcome == "lumBHER2neg_data", "Luminal B/HER2-negative like",
                                               ifelse(all_presso$Outcome == "lumB_data", "Luminal B-like",
                                                      ifelse(all_presso$Outcome == "HER2enr_data", "HER2-enriched-like", "Triple negative")))))
# Check
table(all_presso$Outcome, all_presso$Outcome_def) # correct


# Merge MR PRESSO to table 2 object
all_riskfactors_table2_inclpresso <- merge(x = all_riskfactors_table2,
                                           y = all_presso[, c("Risk_factor", "Outcome_def", "PRESSO.est", "PRESSO.se")],
                                           by.x = c("Risk factor", "Outcome"),
                                           by.y = c("Risk_factor", "Outcome_def"),
                                           all.x = T)


# Combine data for all risk factors into one dataframe
all_riskfactors_fig1 = rbind(height_all_fig1,
                             bmi_all_fig1,
                             t2d_all_fig1,
                             menarche_all_fig1,
                             menopause_all_fig1,
                             density_all_fig1,
                             alcohol_all_fig1,
                             smoking_all_fig1,
                             activity_all_fig1)


# Subset data to results including PCs only
fig1_pcs = subset(all_riskfactors_fig1,
                  all_riskfactors_fig1$Analysis == "IVW including PCs")


# Merge IVW est and se's to input_supplementfig
# First change names of these columns
colnames(fig1_pcs)[4] <- "prim.IVW.est"
colnames(fig1_pcs)[5] <- "prim.se"

all_riskfactors_table2_ib <- merge(x = all_riskfactors_table2_inclpresso,
                                   y = fig1_pcs[, c("Risk factor", "Outcome", "prim.IVW.est", "prim.se")],
                                   by.x = c("Risk factor", "Outcome"),
                                   by.y = c("Risk factor", "Outcome"),
                                   all.x = T)


# For table 2: calculate ORs and 95%CIs for each method + select relevant columns
table2_manuscript <- all_riskfactors_table2_ib %>%
  arrange(factor(`Risk factor`, levels = unique(all_riskfactors_table2$`Risk factor`)), factor(Outcome, levels = (unique(all_riskfactors_table2$Outcome)))) %>%
  select(c(`Risk factor`, Outcome, prim.IVW.est, prim.se, IVW.est, IVW.se, Egger.est, Egger.se, P.pleiotr, Median.est, Median.se, Mode.est, Mode.se, PRESSO.est.y, PRESSO.se.y)) %>%
  mutate(prim.IVW.CI = paste(round(exp(prim.IVW.est - (1.96*prim.se)), digits = 2), round(exp(prim.IVW.est + (1.96*prim.se)), digits = 2), sep = ", ")) %>%
  mutate(IVW.CI = paste(round(exp(IVW.est - (1.96*IVW.se)), digits = 2), round(exp(IVW.est + (1.96*IVW.se)), digits = 2), sep = ", ")) %>%
  mutate(Egger.CI = paste(round(exp(Egger.est - (1.96*Egger.se)), digits = 2), round(exp(Egger.est + (1.96*Egger.se)), digits = 2), sep = ", ")) %>%
  mutate(Median.CI = paste(round(exp(Median.est - (1.96*Median.se)), digits = 2), round(exp(Median.est + (1.96*Median.se)), digits = 2), sep = ", ")) %>%
  mutate(Mode.CI = paste(round(exp(Mode.est - (1.96*Mode.se)), digits = 2), round(exp(Mode.est + (1.96*Mode.se)), digits = 2), sep = ", ")) %>%
  mutate(PRESSO.CI = paste(round(exp(PRESSO.est.y - (1.96*PRESSO.se.y)), digits = 2), round(exp(PRESSO.est.y + (1.96*PRESSO.se.y)), digits = 2), sep = ", ")) %>%
  mutate(across(matches(".est"), exp)) %>%
  mutate(across(matches(".est"), round, 2)) %>%
  select(c(`Risk factor`, Outcome, prim.IVW.est, prim.IVW.CI, IVW.est, IVW.CI, Egger.est, Egger.CI, P.pleiotr, Median.est, Median.CI, Mode.est, Mode.CI, PRESSO.est.y, PRESSO.CI))


table2_manuscript_def <- table2_manuscript %>%
  mutate(across(matches(".CI"), ~ paste0(.x, ")"))) %>%
  mutate(across(matches(".CI"), ~paste0("(", .x)))

print(table2_manuscript_def, quote = T)

# For figures: Subset new table 2 to relevant estimates
all_riskfactors_table2_ib_2 <- all_riskfactors_table2_ib %>%
  select(c(`Risk factor`, Outcome, prim.IVW.est, prim.se, IVW.est, IVW.se, Egger.est, Egger.se, Median.est, Median.se, Mode.est, Mode.se, PRESSO.est.y, PRESSO.se.y))


input_supplementfig <- reshape(data = all_riskfactors_table2_ib_2,
                               idvar = c("Risk factor", "Outcome"),
                               varying = list(c("prim.IVW.est", "IVW.est", "Egger.est", "Median.est", "Mode.est", "PRESSO.est.y"), 
                                              c("prim.se", "IVW.se", "Egger.se", "Median.se", "Mode.se", "PRESSO.se.y")),
                               v.names = c("beta", "se"),
                               times = c("Primary IVW estimates (correlated)", "IVW uncorrelated SNPs", "MR-Egger", "Weighted Median", "Weighted Mode", "MR-PRESSO"),
                               direction = "long")

input_supplementfig$order <- ifelse(input_supplementfig$time == "Primary IVW estimates (correlated)", 1,
                                    ifelse(input_supplementfig$time == "IVW uncorrelated SNPs", 2,
                                           ifelse(input_supplementfig$time == "MR-Egger", 3,
                                                  ifelse(input_supplementfig$time == "Weighted Median", 4,
                                                         ifelse(input_supplementfig$time == "Weighted Mode", 5, 6)))))

rm(all_riskfactors_table2_ib,
   all_riskfactors_table2_ib_2)


risk_factors <- unique(input_supplementfig$`Risk factor`)

plots <- list()

for (i in risk_factors){
  
  plots[[i]] <- input_supplementfig %>%
    filter(`Risk factor` == i) %>%
    ggplot(aes(y = reorder(time, -order), x = exp(beta))) +
    geom_point(size = 1.5) + # check, the size of the point for odd ratio
    coord_trans(x = "log10") +
    scale_x_continuous(limits = c(0.001, 6), breaks = c(0.1, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2, 3, 4, 5, 6)) +
    facet_wrap(~factor(Outcome, levels = (unique(all_riskfactors_table2$Outcome))), nrow = 6, scales = "fixed") +
    theme_bw() +
    theme(legend.position="none") +
    geom_segment(aes(x = exp(beta - (1.96*se)), xend = exp(beta + (1.96*se)), yend = time), size = 0.5) +
    geom_vline(lty = 2, aes(xintercept = 1), colour = "grey") +
    xlab("Causal odds ratio \n\nSupplemental Figure 2H. Forrest plots with causal effects of smoking on each hormone receptor breast cancer subtype across primary and secondary MR methods.") + # needs to be changed for each risk factor
    ylab("") +
    ggtitle(i) +
    theme(axis.text.x = element_text(size = 9), # the numbers of odd ratio
          axis.title.x = element_text(size = 10)) + # the title of odd ratio
    theme(axis.text.y = element_text(size = 9)) + # size of the text in y axis
    theme(strip.text.x = element_text(size = 9)) # label of the box
  
} # note to self need to check estimates with earlier version table 2 (including ORs + CIs)

plot_menarche <- plots[[1]] 
plot_menarche 

plot_menopause <- plots[[2]]
plot_menopause

plot_alcohol <- plots[[3]]
plot_alcohol

plot_bmi <- plots[[4]]
plot_bmi

plot_density <- plots[[5]]
plot_density

plot_height <- plots[[6]]
plot_height

plot_activity <- plots[[7]]
plot_activity

plot_smoking <- plots[[8]]
plot_smoking

plot_t2d <- plots[[9]]
plot_t2d



############################# supplemental figure 3 ############################
# Causal effect estimates for age at menarche from univariable and multivariable MR analysis.

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

# Format menarche estimates in a way that are appropriate for plotting
menarche_ivw <- menarche[, c(1:5)]
menarche_ivw$order <- rep(1:6)
menarche_ivw$analysis <- "Univariable MR"

menarche_mv_yengo <- menarche[, c(1:3, 8,9)]
menarche_mv_yengo$order <- rep(1:6)
menarche_mv_yengo$analysis <- "Multivariable MR, Yengo et al."
names(menarche_mv_yengo) <- names(menarche_ivw)


menarche_mv_ukb <- menarche[, c(1:3, 10,11)]
menarche_mv_ukb$order <- rep(1:6)
menarche_mv_ukb$analysis <- "Multivariable MR, UKB"
names(menarche_mv_ukb) <- names(menarche_ivw)


menarche_plot <- rbind(menarche_ivw,
                       menarche_mv_yengo,
                       menarche_mv_ukb)


# Create plot stratified by subtype
multivar_plot <- ggplot(menarche_plot, aes(y = analysis, x = exp(IVW.est), color = analysis)) +
  geom_point(size = 4) +
  coord_trans(x = "log10") +
  scale_colour_manual("IVW analysis", values = c(inauguration("inauguration_2021_bernie")[c(1, 6, 4)])) +
  #, breaks = c("TN", "LumB", "LumA", "HER2")) +
  facet_wrap(~reorder(Outcome, order), ncol = 6, scales = "fixed") +
  theme_bw() +
  theme(legend.position="none") +
  geom_segment(aes(x = exp(IVW.est - (1.96*IVW.se)), xend = exp(IVW.est + (1.96*IVW.se)), yend = analysis), size = 1) +
  geom_vline(lty = 2, aes(xintercept = 1), colour = "grey") +
  xlab("Causal odds ratio") +
  ylab("") +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold"))

multivar_plot

# Create variable to adjust order when results are stratified by analysis type
menarche_plot$analysis_order <- ifelse(menarche_plot$analysis == "Univariable MR", 1,
                                       ifelse(menarche_plot$analysis == "Multivariable MR, Yengo et al.", 2, 3))

multivar_plot_2 <- ggplot(menarche_plot, aes(y = reorder(Outcome, -order), x = exp(IVW.est), color = Outcome)) +
  geom_point(size = 4) +
  coord_trans(x = "log10") +
  scale_colour_manual("Molecular subtype", values = c(inauguration("inauguration_2021_bernie")[1:6])) +
  #, breaks = c("TN", "LumB", "LumA", "HER2")) +
  facet_wrap(~reorder(analysis, analysis_order), ncol = 3, scales = "fixed") +
  theme_bw() +
  theme(legend.position="none") +
  geom_segment(aes(x = exp(IVW.est - (1.96*IVW.se)), xend = exp(IVW.est + (1.96*IVW.se)), yend = Outcome), size = 1) +
  geom_vline(lty = 2, aes(xintercept = 1), colour = "grey") +
  xlab("Causal odds ratio") +
  ylab("") +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold"))

multivar_plot_2


## the end of the script ##

