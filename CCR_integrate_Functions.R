# This is code contains functions to integrate CCR copy number and gene expression data.
# Code was developed by Samir Lal 



###function to merge expression data files normalise quantiles and remove batch effects COMBAT 
get.gene.rna<-function(GEX1,GEX2,probes.gex1,probes.gex2,genes){
  colnames(GEX1)<-gsub(".AVG_Signal","",colnames(GEX1))
  colnames(GEX2)<-gsub(".AVG_Signal","",colnames(GEX2))

  GEX1<-GEX1[,colnames(GEX1) %in% ILC$array_name]
  GEX2<-GEX2[,colnames(GEX2) %in% ILC$array_name]

  colnames(GEX1)<-ILC[match(colnames(GEX1),ILC$array_name),1]
  colnames(GEX2)<-ILC[match(colnames(GEX2),ILC$array_name),1]



  rownames(GEX1)<-probes.gex1
  rownames(GEX2)<-probes.gex2

  ####log transform the raw expression values
  
  eset1<-data.matrix(log2(GEX1))
  eset2<-data.matrix(log2(GEX2))

  ### remove Poor quality samples
  eset2<-eset2[,!(colnames(eset2) %in% c("L18"))]

  #Apply quantile normalisation
  eset1<-normalizeQuantiles(eset1)
  eset2<-normalizeQuantiles(eset2)
  pdat1<-ILC[match(colnames(eset1),ILC$ID),]
  pdat1$Study<-rep(1,nrow(pdat1))
  pdat2<-ILC[match(colnames(eset2),ILC$ID),]
  pdat2$Study<-rep(2,nrow(pdat2))
  
  row.names(pdat1)<-pdat1$ID
  row.names(pdat2)<-pdat2$ID
  
  #Create expression data set
  eset1 <- ExpressionSet(as.matrix(eset1),phenoData = AnnotatedDataFrame(data=pdat1))
  eset2 <- ExpressionSet(as.matrix(eset2),phenoData = AnnotatedDataFrame(data=pdat2))

  esets<-list(eset1,eset2)


  #merge the expression datasets without batch effect correction 
  eset_FRMA<-inSilicoMerging::merge(esets)
  #Plot multidimensional scaling plot (Distance between expression profiles)
  plotMDS(eset_FRMA,colLabel="Hist.Diagnosis",symLabel="Study",main="FRMA (No Transformation)",file="/Users/uqslal/Desktop/our_data_batch_MDS.pdf")
  #select the samples
  select = sample(1:ncol(eset_FRMA),25);
  #Plot the relative log expression
  plotRLE(eset_FRMA[,select], colLabel="Study", legend=FALSE, main="FRMA",file="/Users/uqslal/Desktop/our_data_batch_RLE.pdf");

  #merge the expression datasets with batch effect correction  
  eset_COMBAT<-inSilicoMerging::merge(esets,method="COMBAT")
  #Plot multidimensional scaling plot (Distance between expression profiles)
  plotMDS(eset_COMBAT,colLabel="Hist.Diagnosis",symLabel="Study",main="COMBAT",file="/Users/uqslal/Desktop/COMBAT_our_data_batch_MDS.pdf")
  #select the samples
  select = sample(1:ncol(eset_COMBAT),25);
  #Plot the relative log expression 
  plotRLE(eset_COMBAT[,select], colLabel="Study", legend=FALSE, main="COMBAT",file="/Users/uqslal/Desktop/COMBAT_our_data_batch_RLE.pdf");


  exprs<-as.data.frame(exprs(eset_COMBAT))

  #GEX1<-read.table("/Volumes/DISK1/SampleProbeProfile_200912.txt",sep="\t",header=TRUE)

  exprs$genes<-genes[match(row.names(exprs),genes$PROBE_ID),2]
return(exprs)
}



#function to read in the copy number segmented data and the 
genes.copy.number<-function(tss,UQCCR,logR=FALSE){
  gr0 = with(tss, GRanges(V3, IRanges(start=V5, end=V6)))
  gr1 = with(UQCCR, GRanges(chrom, IRanges(start=start.pos, end=end.pos)))
  hits = findOverlaps(gr0, gr1,type="within")
  
  res<-cbind(UQCCR[subjectHits(hits),],tss[queryHits(hits),5])
  
  colnames(res)[ncol(res)]<-"genes"
  colnames(res)[7]<-'CN'
#Change the copy number matrix from to numeric states 1: gain,2:amp,0: neut,-1:hetd,-2:homd
  res$V8<-res$CN
  res$CN<-NULL
  neg<-function(x)-x
  res[res$V8>3 & res$V8<=5,9]<-1
  res[res$V8>5,9]<-2
  res[res$V8<2 & res$V8>0 ,9]<-neg(1)
  res[res$V8<1,9]<-neg(2)
  res[is.na(res$V9),9]<-0
  res$V8<-NULL
  res$sampleID<-as.character(res$sampleID)
  res.list<-split(res,f=res$sampleID)
  
  for(i in 1:length(res.list)){
    res.list[[i]]<-res.list[[i]][order(abs(res.list[[i]][,8]),abs(res.list[[i]][,6]),decreasing=TRUE),]
    res.list[[i]]<-res.list[[i]][!duplicated(res.list[[i]][,7]),]
  }
  
  mydf <- do.call("rbind",res.list)
  #### 2 options one will generate gene matrix with LogR values (for spearman correlation) 
  #### The other will generate a gene matrix with raw copy number values (for ANOVA)      
  if(logR==FALSE){
    mydf<-mydf[,c(1,2,3,4,5,8,7)];
  }
  if(logR==TRUE){
    mydf<-mydf[,c(1,2,3,4,5,6,7)] 
  }
  
  mydf.list<-split(mydf[,c(6,7)],f=mydf$sampleID)
  
  
  for(i in 1:length(mydf.list)){
    colnames(mydf.list[[i]])[1]<-names(mydf.list)[i]
  }
  #library(plyr); is used here as well as several inbuilt functions i.e. join_all.
  out<-join_all(mydf.list,by="genes")
  
  rownames(out)<-out$genes
  out$genes<-NULL
  return(out)
}



####merge GEX and CN matrix priming for integration
combine_CN_GEX<-function(gene_mat,rna_mat){
  genes<-rna_mat$genes;
  rna_mat<-rna_mat[,which(colnames(rna_mat) %in% colnames(gene_mat))];
  gene_mat<-gene_mat[,colnames(gene_mat) %in% colnames(rna_mat)];
  rna_mat<-rna_mat[,colnames(gene_mat)];
  rna_mat$genes<-genes
  eset<-base::merge(gene_mat,rna_mat,by.x="row.names",by.y="genes");
  return(eset);
}



integrate_anova<-function(ILC_matrix){
  neg<-function(x)-x
  CN<-apply(ILC_matrix[,c(2:((ncol(ILC_matrix)-1)/2)+1)],2, function(x) ifelse(as.numeric(as.character(x))>1.1,"AMP",x))
  CN<-apply(CN,2, function(x) ifelse(as.numeric(as.character(x))>0.5 & x!="AMP", "GAIN",x))
  CN<-apply(CN,2, function(x) ifelse(as.numeric(as.character(x))<neg(1.1) & x!="AMP" & x !="GAIN", "HOMD",x))
  CN<-apply(CN,2, function(x) ifelse(as.numeric(as.character(x))<neg(0.5) & x!="AMP" & x !="GAIN" & x!="HOMD", "LOSS",x))
  CN<-apply(CN,2, function(x) ifelse(x=="0","NEUT",x))
  genes<-ILC_matrix[,1]
  GEX<-ILC_matrix[,c((((ncol(ILC_matrix)-1)/2)+2):ncol(ILC_matrix))];
  #####Apply the anova
  #vector for the p-values
  pvals<-c();
  #vector for the change that is contributing to the p-value 
  changes<-c();
  #vector of YES or NO for significant cases.
  significant<-c()
  #vector for the number of samples
  samples<-c()
  #vector with the variances will be used in the meta-analysis
  vars<-c()
  for (i in 1:nrow(GEX)){
    out<-t(rbind(GEX[i,],CN[i,]));
    dats<-as.data.frame(out);
    colnames(dats)<-c("GEX","CN");
    dats<-dats[!is.na(dats$CN),];
    dats<-dats[!is.na(dats$GEX),];
    samples<-c(samples,nrow(dats));
    dats$GEX<-as.numeric(as.character(dats$GEX));
    dats$GEX<-as.numeric(as.character((dats$GEX)));
    if(length(unique(dats$CN))>1){
      fit<-lm(GEX~CN,data=dats);
      anova.res<-anova(fit);
      ao1<-aov(GEX~CN,data=dats);
      #Post hoc test to determine where significant change lies
      result.t<-TukeyHSD(ao1,"CN");
      change<-result.t$CN;
      changes<-c(changes,row.names(change[order(-change[,4]), , drop = FALSE])[1]);
      ##Is my result significant P-value <0.05
      sig<-result.t$CN[,4][which(result.t$CN[,4]<0.05)];
      sig<-ifelse(length(sig)<1,"FALSE","TRUE");
      significant<-c(significant,sig);
      pval<-anova.res$"Pr(>F)"[1];
      pvals<-c(pvals,pval);
    }
    else{
      pval<-NA;
      significant<-c(significant,NA);
      changes<-c(changes,NA);
      pvals<-c(pvals,pval);
    }
    
    
  }
  ILC_matrix<-cbind(ILC_matrix,pvals,significant,changes,samples)
  return(ILC_matrix)
}




####Function to integrate using Spearman 
integrate_correlation<-function(ILC_matrix){
  #store the number of samples
  samples.out<-apply(ILC_matrix[,2:(ncol(ILC_matrix)/2)],1,function(x) length(which(!is.na(x))));
  #store the rho
  rho.out<-apply(ILC_matrix[,(2:ncol(ILC_matrix))],1,function(x) return(cor.test(x[c(1:((length(x)-1)/2))],x[c((((length(x)-1)/2)+2):length(x))],method="spearman",alternative="greater",na.action=FALSE)$estimate))
  #store the p-values
  pval.out<-apply(ILC_matrix[,(2:ncol(ILC_matrix))],1,function(x) return(cor.test(x[c(1:((length(x)-1)/2))],x[c((((length(x)-1)/2)+2):length(x))],method="spearman",alternative="greater",na.action=FALSE)$p.val))
  #combine and collate this data
  out<-as.data.frame(cbind(ILC_matrix$Row.names,rho.out,pval.out,samples.out))
  return(out)
}
