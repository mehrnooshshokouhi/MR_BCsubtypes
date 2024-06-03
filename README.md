# Causal effects of breast cancer risk factors across hormone receptor breast cancer subtypes: A two-sample Mendelian randomization study
This repository includes the codes used to perform analyses for the mentioned study.   
# Part 1: Rscripts

####1. MR_LDlink.R###
   
   The aim of this Rscript is to identify BCAC proxies for MR analyses.
   
   For SNPs that were not available in the BCAC GWAS summary statistics, we searched for proxy SNPs (linkage disequilibrium r2 ≥ 0.8) using the NIH LDlink API implemented in the LDlinkR R-package.
   
   To search for proxy SNPs for each risk factor, we used the following functions:
   
   1.1 Function to make variables to indicate index SNPs and overlap with BCAC data

   1.2 Function to select one BCAC proxy for each index SNP based on the highest R2

   1.3 Function to select one BCAC proxy for each index SNP based on the highest R2, excluding indels;
       Only use this function if none of the index variants are indels!

###2. MR_Objects_Riskfactors.R###
   
   The aim of this Rscript is to make final risk factor data frames for MR analysis.

   Prior to conducting MR analyses, we performed harmonization of alleles and effect estimates between the risk factor and BCAC GWAS using the TwoSampleMR R-package.

   At this step, we excluded palindromic SNPs with intermediate allele frequencies (i.e., A/T or C/G SNPs with an effect allele frequency ranging from 0.40 to 0.60) because harmonization of these specific variants     between different data sources is very error prone.

   In addition, SNPs that were not available in the 1000G phase 3 reference panel were excluded during harmonization.

   To finalized the genetic instrumental variables, we used the following functions:
   
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
   
3. MR_Primary_And_Secondary_Analysis.R
   
6. MR_PRESSO_Analysis.R
   
7. MR_Manuscript_And_Supplemental_Figures.R
    
8. MR_Manuscript_And_Suuplemental_Tables.R
    
# Part 2: Summary-level data
