LobSig Scripts that perform integrative analysis and LobSig scoring.

Make sure you change the directories of the output files for your local machine.  

Large data files .RData  Expression_CopyNumber_files.RData contact the corresponding authors.

0-master_Rscript is the only file that needs to be executed this file loads 4 files and performs all the neccessary tests ANOVA and Spearman with gene expression and copy number datasets 
CCR_integrate_Functions.R which contains the functions necessary to integrate CCR dataset 
TCGA_integrate_Functions.R which contains the functions necessary to integrate TCGA dataset 
METABRIC_integrate_Functions.R which contains the functions necessary to integrate METABRIC dataset 
combine_p_values.R which contains the functions necessary for the meta-analysis 

MachineLearning.R contains the functions necessary to normalize gene expression data and provide the LobSig score and classifications as "High" or Low". 
