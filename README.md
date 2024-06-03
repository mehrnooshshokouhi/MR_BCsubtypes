# Causal effects of breast cancer risk factors across hormone receptor breast cancer subtypes: A two-sample Mendelian randomization study
This repository includes the codes used to perform analyses for the mentioned study.   
## Part 1: Rscripts

### 1. MR_LDlink.R 
   
   * The aim of this Rscript is to identify BCAC proxies for MR analyses.
   
      * For SNPs that were not available in the BCAC GWAS summary statistics, we searched for proxy SNPs (linkage disequilibrium r2 ≥ 0.8) using the NIH LDlink API implemented in the LDlinkR R-package.
   
      #### To search for proxy SNPs for each risk factor, we used the following functions:
   
      1.1 Function to make variables to indicate index SNPs and overlap with BCAC data

      1.2 Function to select one BCAC proxy for each index SNP based on the highest R2

      1.3 Function to select one BCAC proxy for each index SNP based on the highest R2, excluding indels; 
          Only use this function if none of the index variants are indels!

### 2. MR_Objects_Riskfactors.R
   
   * The aim of this Rscript is to make final risk factor data frames for MR analysis.

      * Prior to conducting MR analyses, we performed harmonization of alleles and effect estimates between the risk factor and BCAC GWAS using the TwoSampleMR R-package.

      * At this step, we excluded palindromic SNPs with intermediate allele frequencies (i.e., A/T or C/G SNPs with an effect allele frequency ranging from 0.40 to 0.60) because harmonization of these specific              variants between different data sources is very error prone.

      * In addition, SNPs that were not available in the 1000G phase 3 reference panel were excluded during harmonization.

      #### To finalized the genetic instrumental variables, we used the following functions:
   
      2.1 Function for creating indicator variable proxy_needed

      2.2 Function to separate correlated alleles for index and proxy SNPs; using a two-step approach
  
      2.3 Function to make indicator variable for palindromic SNPs (index SNPs)
  
      2.4 Function to make indicator variable for palindromic SNPs (proxy SNPs)
  
      2.5 Function to match alleles
  
      2.6 Function to check if matching of alleles went correctly
   
      2.7 Function to make final snp id variable

      2.8 Function to include reason for absence proxy SNP

      2.9 Function to create supplemental table 1

      2.10 Function to create final MR risk factor object
   
### 3. MR_Primary_And_Secondary_Analysis.R

   * This Rscript includes the final two-sample Mendelian randomization analysis for the manuscript.
     
      * Prior to performing our primary MR analyses, we set out to calculate linkage disequilibrium matrices for all SNPs for the specific risk factors using the ld_matrix function implemented in the TwoSampleMR           R-package.

      * However, due to the substantial proportion of highly correlated SNPs for height and BMI, the correlation matrices that we calculated for these risk factors were near-singular. Therefore, we performed               unscaled principal component analysis on a weighted version of the genetic correlation matrix instead. 

      * We employed the IVW method using a multiplicative random-effects model to calculate primary MR estimates for all nine risk factors in relation to each hormone receptor breast cancer subtype. 

      * We additionally performed IVW analyses restricted to uncorrelated SNPs (linkage disequilibrium r2 ≤ 0.001).

      #### To perform the primary and secondary MR analysis, we used the following functions:

      ##### 3.1 Create a first MR pipeline that:
     
      3.1.1 Harmonizes the exposure-outcome data
         
      3.1.2 Converts the harmonized data to an MRinput object (includes calculation LD matrix)
         
      3.1.3 Performs a IVW analysis with multiplicative effects including a LD matrix
         
      3.1.4 Saves the results from these analyses

      ##### 3.2 Create a third MR pipeline that:
     
      3.2.1 Clumps the genetic instruments
         
      3.2.2 Harmonizes the clumped exposure-outcome data
         
      3.2.3 Converts the harmonized data to an MRinput object
         
      3.2.4 Performs an IVW analysis with multiplicative effects
         
      3.2.5 Performs a MR-Egger analysis
 
      3.2.6 Performs a weighted median analysis
         
      3.2.7 Performs a weighted mode analysis
         
      3.2.8 Performs a MR-PRESSO analysis
         
      3.2.9 Saves the results from these analyses
   
### 4. MR_PRESSO_Analysis.R

   * The aim of this Rscript is to perform the MR-PRESSO analyses for all nine risk factors.
     
      * MR_PRESSO was performed via R base in the terminal. Because of this, it has a separate script.
   
### 5. MR_Manuscript_And_Supplemental_Figures.R

   * The aim of this Rscript is to create figures in the main manuscript, and in the supplemental file, including:

     manuscript figure 1. Causal breast cancer subtype-specific effect estimates per unit increase for nine established breast cancer risk factors. 
    
     manuscript figure 2. Overview of evidence for (subtype-specific) causal effects per increasing unit of the risk factor. 
    
     supplemental figure 2. Forrest plots with causal effects of nine risk factors on each hormone receptor breast cancer subtype across primary and secondary MR methods.
    
     supplemental figure 3. Causal effect estimates for age at menarche from univariable and multivariable MR analysis.
 
### 6. MR_Manuscript_And_Suuplemental_Tables.R
    
## Part 2: Summary-level data
