# This is code contains functions to integrate METABRIC copy number and gene expression data.
# Code was developed by Samir Lal 




####Function to read in the segmented copy number files and refseq gene file  
####Generates a copy number x sample matrix 
get.genes.copynumber<-function(segmented_data,tss,logR=FALSE){
  tss<-tss[,c(3,5,6,2,13)];
  tss$V3<-gsub("chrX","chr23",tss$V3);
  tss$V3<-gsub("chr","",tss$V3);
  print("Loading Genomic Ranges.....");
  #using library(GenomicRanges) and the inbuilt functions;
  gr0 = with(tss, GRanges(V3, IRanges(start=V5, end=V6)));
  gr1 = with(segmented_data, GRanges(chrom, IRanges(start=start.pos, end=end.pos)));
  print("overlapping");
  hits = findOverlaps(gr0, gr1,type="within")
  res<-cbind(segmented_data[subjectHits(hits),],tss[queryHits(hits),5]);
  colnames(res)[ncol(res)]<-"genes";
  res$V8<-res$CN
  res$CN<-NULL
  res$sampleID<-as.character(res$sampleID)
  res.list<-split(res,f=res$sampleID)
  for(i in 1:length(res.list)){
    res.list[[i]]<-res.list[[i]][order(abs(res.list[[i]][,8]),abs(res.list[[i]][,6]),decreasing=TRUE),]
    res.list[[i]]<-res.list[[i]][!duplicated(res.list[[i]][,7]),]
    sample<-i;
    print(paste(paste("Processing sample",sample,sep=""),paste("of",length(res.list),sep=""),sep=""))
  }
  print("assembling processed samples");
  mydf <- do.call("rbind",res.list);
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
  print("Loading plyr.....");
  #library(plyr); is used here as well as several inbuilt functions i.e. join_all.
  out<-join_all(mydf.list,by="genes")
  rownames(out)<-out$genes
  out$genes<-NULL
  return(out)
  }



#####Function to read in the raw RNA-seq expression files and generate a gene x sample matrix of expression values

get.tumour.genes.RNA<-function(exprs.disc,exprs.valid,out){
  exprs<-base::merge(exprs.disc,exprs.valid,by="row.names");
  row.names(exprs)<-exprs$Row.names;
  exprs$Row.names<-NULL;
  genes<-exprs$GENE.SYMBOL;
  ILC<-exprs[,colnames(exprs) %in% colnames(out)];
  ILC<-ILC[,colnames(out)];
  #quantile normalize discovery and validation ILC datasets. 
  ILC<-normalizeQuantiles(ILC);
  ILC$genes<-genes;
  return(ILC)
}


####merge GEX and CN matrix priming for integration
merge_CN_GEX<-function(CN,GEX){
  ILC<-base::merge(CN,GEX,by.x="row.names",by.y="genes");
  return(ILC)
}


####Function to integrate using ANOVA
integrate_anova<-function(ILC_matrix){
  ####Change the copy number matrix from numeric states to characters
  neg<-function(x)-x
  CN<-apply(ILC_matrix[,c(2:((ncol(ILC_matrix)-1)/2))],2, function(x) ifelse(as.numeric(as.character(x))>1.1,"AMP",x))
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




