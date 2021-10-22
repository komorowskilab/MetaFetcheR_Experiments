# MetaFetcheR_Experiments
This repo includes scripts for the reproducing figures used in the manuscript and the comparisons:

The datasets are in the 2 excel files :
1-20190426_metabolite_ids.xlsx (Diamanti et al)
2-prostate_cancer_dataset.xlsx (Priolo et al)

The data used to generate figure 4 can be loaded using 
ResultsCoverage.RData to generate the figure 

You will need to use the function in "utils.R" to run the scripts 

You need to install and load the required packages

library(MetaboAnalystR)
library(metafetcher)
library(xlsx)
library(dplyr)
library(ggplot2)
library(ggsci)
library(wesanderson)


The scripts used for the comparisons of MetaFetcheR with MS_Targeted and MetaboAnalystR on both datasets (Diamanti et al and Priolo et al) (case 1 and case2)  can be found in the script Comparisons_MetaFetcheR_MetaboAnalystR_MS_Targeted.R

The scripts used for the comparison of MetaFetcheR with CTS (case 3)can be found in the script Comparison_MetaFetcheR_CTS.R
The source of data is in  CTS_comparisons/CTSInputMetaFetcheR.xlsx

The script used for the comparison of MetaFetcheR with CTS (case 3)can be found in the script  Comparison_MetaFetcheR_MA5.R 
The source of data is in  MA5Comparison/MA5InputMetaFetcheR.xlsx

The source of the data for the figures can be found in Source1.xlsx and Source2.xlsx 


Comments are included in each section of all the scripts to guide the user

