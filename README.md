LobSig Scripts that perform integrative analysis and LobSig scoring as featured in the LobSig manuscript.

Make sure you change the directories of the output files for your local machine.  

Contact the corresponding authors for large data files .RData  Expression_CopyNumber_files.RData .

The github repository contains the following key scripts.

Required R packages include
genefu https://bioconductor.org/packages/release/bioc/html/genefu.html
caret https://topepo.github.io/caret/
matrixStats https://cran.rstudio.com/web/packages/matrixStats/index.html
ROCR https://rocr.bioinf.mpi-sb.mpg.de/
openxlsx https://cran.r-project.org/web/packages/openxlsx/
GenomicRanges https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
plyr https://www.rdocumentation.org/packages/plyr/versions/1.8.4
limma https://bioconductor.org/packages/release/bioc/html/limma.html
stringr https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html
edgeR https://bioconductor.org/packages/release/bioc/html/edgeR.html
inSilicoMerging https://www.bioconductor.org/packages//2.10/bioc/html/inSilicoMerging.html
Biobase http://bioconductor.org/packages/release/bioc/html/Biobase.html
biomaRt https://bioconductor.org/packages/release/bioc/html/biomaRt.html

0-master_Rscript is the only file that needs to be executed this file loads 4 files and performs all the neccessary tests A
ANOVA and Spearman with gene expression and copy number datasets 

CCR_integrate_Functions.R which contains the functions necessary to integrate CCR dataset 

TCGA_integrate_Functions.R which contains the functions necessary to integrate TCGA dataset 

METABRIC_integrate_Functions.R which contains the functions necessary to integrate METABRIC dataset 
combine_p_values.R which contains the functions necessary for the meta-analysis 

MachineLearning.R contains the functions necessary to normalize gene expression data and provide the LobSig score and classifications as "High" or Low". 
