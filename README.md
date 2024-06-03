# Causal effects of breast cancer risk factors across hormone receptor breast cancer subtypes: A two-sample Mendelian randomization study
This repository includes the codes used to perform analyses for the mentioned study.   
# Part 1: Rscripts
1. MR_LDlink.R
   
   The aim of this Rscript is to identify BCAC proxies for MR analyses.
   
   For SNPs that were not available in the BCAC GWAS summary statistics, we searched for proxy SNPs (linkage disequilibrium r2 ≥ 0.8) using the NIH LDlink API implemented in the LDlinkR R-package.
   
   To search for proxy SNPs for each risk factor, we used three functions:
   
   1.1 Function to make variables to indicate index SNPs and overlap with BCAC data

   1.2 Function to select one BCAC proxy for each index SNP based on the highest R2

   1.3 Function to select one BCAC proxy for each index SNP based on the highest R2, excluding indels;
       Only use this function if none of the index variants are indels!
   
3. MR_Objects_Riskfactors.R
   
4. MR_Primary_And_Secondary_Analysis.R
   
5. MR_PRESSO_Analysis.R
   
6. MR_Manuscript_And_Supplemental_Figures.R
    
7. MR_Manuscript_And_Suuplemental_Tables.R
    
# Part 2: Summary-level data
