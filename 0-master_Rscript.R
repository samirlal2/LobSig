# This is code to replicate the integrative analysis from 2019 NPJ
# Breast cancer paper. Code was developed by Samir Lal

####load required packages 
require(GenomicRanges);
require(plyr);
require(limma);
require(stringr)
require(edgeR);
require(inSilicoMerging)
require(Biobase)
require(biomaRt)


# 3 R scripts with appropriate functions are provided. 
# METABRIC_integrate_Functions.R
# CCR_integrate_Functions.R
# TCGA_integrate_Functions.R

# Set the working directory and output directory. 
wkdir <- "/Users/samirlal/Desktop/renpjbreastcancernpjbcancer00372majorrevision/"
setwd(wkdir)
load("Expression_CopyNumber_files.RData")


  
#####Analyse METABRIC ILC dataset################
##                                             ##
#################################################
source("METABRIC_integrate_Functions.R")
  
  gene_matrix<-get.genes.copynumber(METABRIC,tss.METABRIC,logR=FALSE)
  
  
  gene_matrix1<-get.genes.copynumber(METABRIC,tss.METABRIC,logR=TRUE)
  
  
  RNA_matrix<-get.tumour.genes.RNA(exprs.disc,exprs.valid,gene_matrix)
  
  ILC.merged<-merge_CN_GEX(gene_matrix,RNA_matrix)
  
  ILC.spearman<-merge_CN_GEX(gene_matrix1,RNA_matrix)
  
  print("This takes 8-10 minutes to run so be patient...grab a coffee")
  mat1<-integrate_anova(ILC.merged)
  
  mat<-integrate_correlation(ILC.spearman)
  mat<-mat[order(as.numeric(as.character(mat$rho.out)),decreasing=TRUE),]
  
  write.table(mat1[,c(1,(ncol(mat1)-3):ncol(mat1))],"METABRIC_anova.txt",sep="\t",col.names=T,row.names=F,quote=F)
  
  write.table(mat,"METABRIC_spearman.txt",sep="\t",col.names=T,row.names=F,quote=F)
  
  #####Analyse TCGA ILC dataset###################
  ##                                             ##
  #################################################
  source("TCGA_integrate_Functions.R")
  
  
   gene_matrix<-get.genes.copynumber(TCGA,tss.TCGA,logR=FALSE)
  
  
  gene_matrix1<-get.genes.copynumber(TCGA,tss.TCGA,logR=TRUE)
  
  
  voom_RNA<-voom.trans(RNA_matrix.tumour,RNA_matrix.normal)

  ILC.TCGA<-merge_CN_GEX(gene_matrix,voom_RNA)
  
  ILC.spearman<-merge_CN_GEX(gene_matrix1,voom_RNA)

  print("This takes 8-10 minutes to run so be patient...grab a coffee")
  mat1<-integrate_anova(ILC.TCGA)
  
  
  
  mat<-integrate_correlation(ILC.spearman)
  mat<-mat[order(as.numeric(as.character(mat$rho.out)),decreasing=TRUE),]
  
  write.table(mat1[,c((ncol(mat1)-3):ncol(mat1))],"TCGA_anova.txt",sep="\t",col.names=T,row.names=T,quote=F)
  
  write.table(mat,"TCGA_spearman.txt",sep="\t",col.names=T,row.names=T,quote=F)
  
  #####Analyse CCR ILC dataset#####################
  ##                                             ##
  ################################################# 
  source("CCR_integrate_Functions.R")
  
    
  probes.gex1<-GEX1$PROBE_ID
  probes.gex2<-GEX2$PROBE_ID
  
  GEX1<-GEX1[,grep(".AVG_Signal",colnames(GEX1))]
  GEX2<-GEX2[,grep(".AVG_Signal",colnames(GEX2))]
  
  
  colnames(ILC)
  
  ILC$array_name<-paste(ILC$"Matrix.(AMR/JJ)",ILC$"Position",sep="_")
  
  ILC<-ILC[,c(1,3,ncol(ILC))]
  
  ILC$array_name<-paste("X",ILC$array_name,sep="")
  
  ILC<-ILC[ILC$Hist.Diagnosis=="ILC",]
  
  rna_matrix<-get.gene.rna(GEX1,GEX2,probes.gex1,probes.gex2,genes)
  
  
  
  ##########CN analysis
  
  tss.CCR<-tss.CCR[,c(3,5,6,2,13)]
  tss.CCR$V3<-gsub("chrX","chr23",tss.CCR$V3)
  
  tss.CCR$V3<-gsub("chr","",tss.CCR$V3)
  
  gene_matrix<-genes.copy.number(tss.CCR,CCR,logR=FALSE)
  
  gene_matrix1<-genes.copy.number(tss.CCR,CCR,logR=TRUE)
  
  ILC<-combine_CN_GEX(gene_matrix,rna_matrix)
  
  ILC.spearman<-combine_CN_GEX(gene_matrix1,rna_matrix)
  
  print("This takes 8-10 minutes to run so be patient...grab a coffee")
  mat1<-integrate_anova(ILC)
  
  mat<-integrate_correlation(ILC.spearman)
  mat<-mat[order(as.numeric(as.character(mat$rho.out)),decreasing=TRUE),]
  
  write.table(mat1[,c(1,(ncol(mat1)-3):ncol(mat1))],"CCR_anova.txt",sep="\t",col.names=T,row.names=T,quote=F)
  
  write.table(mat,"CCR_spearman.txt",sep="\t",col.names=T,row.names=F,quote=F)
  
  ####Perform the meta-analysis across all 3 datasets######
  ##                                                     ##
  #########################################################
  
  source("combine_p-values.R")
  
  ####read in ANOVA results 
  stud1<-read.table("CCR_anova.txt",sep="\t",header=TRUE)  
  stud2<-read.table("METABRIC_anova.txt",sep="\t",header=TRUE) 
  stud3<-read.table("TCGA_anova.txt",sep="\t",header=TRUE,row.names=1)  
  
  #implements Stouffers Z-score
  p.vals.anova<-combine.p.values.anova(stud1,stud2,stud3)
  
  p.vals.anova<-p.vals.anova[order(p.vals.anova$combined.p),]
  
  ###Lets annotate the genes with genomic regions
  
  genes<-p.vals.anova 
  
  g1 <- genes$Row.names
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  # Only use standard human chromosomes
  normal.chroms <- c(1:22, "X", "Y", "M")
  # 
  my.symbols <- g1
  # # Filter on HGNC symbol and chromosome, retrieve genomic location and band
  my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                      filters = c("hgnc_symbol", "chromosome_name"),values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),mart = ensembl)  
  
  Genes<-base::merge(genes,my.regions,by.x="Row.names",by.y="hgnc_symbol")
  
  
  p.vals.anova<-Genes[order(Genes$combined.p),]
  
  gsub(".x",".CCR",colnames(p.vals.anova))
  gsub(".y",".METABRIC",colnames(p.vals.anova))
  
  write.table(p.vals.anova,"anova_combined.txt",sep="\t",col.names=T,row.names=F,quote=F)
  
  stud1<-read.table("CCR_spearman.txt",sep="\t",header=TRUE)  
  stud2<-read.table("METABRIC_spearman.txt",sep="\t",header=TRUE) 
  stud3<-read.table("TCGA_spearman.txt",sep="\t",header=TRUE,row.names=1) 
  
  #implements a random effects model to combine correlation coefficients 
  p.vals.spearman<-combine.p.values.spearman(stud1,stud2,stud3)
  
  p.vals.spearman<-p.vals.spearman[!is.na(p.vals.spearman$efsizes),];
  p.vals.spearman<-p.vals.spearman[order(as.numeric(p.vals.spearman$efsizes),decreasing=TRUE),];
  
  ###Lets annotate the genes with genomic regions 
  
  genes<-p.vals.spearman 
  g1<-row.names(genes)
  

  # # Filter on HGNC symbol and chromosome, retrieve genomic location and band
  my.symbols <- g1
  my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                      filters = c("hgnc_symbol", "chromosome_name"),values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),mart = ensembl)  
  
  p.vals.spearman<-base::merge(genes,my.regions,by.x="row.names",by.y="hgnc_symbol")
  
  
  colnames(p.vals.spearman)<-gsub(".x","CCR",colnames(p.vals.spearman));
  colnames(p.vals.spearman)<-gsub(".y","METABRIC",colnames(p.vals.spearman));
  
  
  
  write.table(p.vals.spearman,"spearman_results.txt",sep="\t",row.names=T,col.names=T,quote=F)
  
  
   
  