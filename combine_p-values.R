

####load required packages
require(biomaRt)
require(meta)


#####Function to combine p-values. 
'combine.test' <-
  function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {
    if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }
    method <- match.arg(method)
    na.ix <- is.na(p)
    if(any(na.ix) && !na.rm) { stop("missing values are present!") }
    if(all(na.ix)) { return(NA) } ## all p-values are missing
    p <- p[!na.ix]
    k <- length(p)
    if(k == 1) { return(p) }
    if(missing(weight)) { weight <- rep(1, k); }
    switch(method,  
           "fisher"={
             cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
           }, 
           "z.transform"={
             z <- qnorm(p, lower.tail=FALSE)
             cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
           }, 
           "logit"={
             tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
             cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
           })
    return(cp)
  }


combine.p.values.anova<-function(stud1,stud2,stud3){
  stud1<-stud1[order(stud1$pvals),];
  stud1<-stud1[!duplicated(stud1$Row.names),];
  #adjust p-values using BH method
  stud1$pvals<-p.adjust(stud1$pvals,method="BH");
  stud2<-stud2[order(stud2$pvals),];
  stud2$pvals<-p.adjust(stud2$pvals,method="BH");
  stud2<-stud2[!duplicated(stud2$Row.names),];
  stud3<-stud3[order(stud3$pvals),];
  stud3$pvals<-p.adjust(stud3$pvals,method="BH");
  stud3<-stud3[!duplicated(row.names(stud3)),];
  merged1<-base::merge(stud1,stud2,by.x="Row.names",by.y="Row.names");
  merged2<-base::merge(merged1,stud3,by.x="Row.names",by.y="row.names");
  zs<-c();
  combine.ps<-c();
  merged2<-merged2[merged2$pvals.x!="NA",];
  merged2<-merged2[merged2$pvals.y!="NA",];
  merged2<-merged2[merged2$pvals!="NA",];
  merged2<-merged2[!is.na(merged2$Row.names),];
   for(i in 1:nrow(merged2)){
    pvals<-merged2[i,c(2,6,10)];
    w<-merged2[i,c(5,9,13)];
    pvals<-as.numeric(pvals);
    ws<-w;
    #implement method that combines p-values
    z<-combine.test(pvals,ws,method="z.transform");
    p<-z[2];
    z<-z[1];
    combine.ps<-c(combine.ps,z);
   }
  merged2$combined.p<-combine.ps
  return(merged2)
}


combine.p.values.spearman<-function(stud1,stud2,stud3){
  stud1<-stud1[!duplicated(stud1$V1),];
  stud2<-stud2[!duplicated(stud2$V1),];
  out<-base::merge(stud1,stud2,by="V1");
  out<-base::merge(out,stud3,by.x="V1",by.y="row.names")
  row.names(out)<-out$V1
  out$V1<-NULL
  cpvals<-c();
  efsizes<-c();
  print(head(out));
  for(i in 1:nrow(out)){
  rhos<-as.numeric(out[i,c(1,4,8)]);
  samples<-as.numeric(out[i,c(3,6,10)]);
  labels<-as.vector(c("CCR","METABRIC","TCGA"));
  if (any(is.na(rhos))){
  cpval<-"NA";
  efs<-"NA";
  cpvals<-c(cpvals,cpval);
  efsizes<-c(efsizes,efs);
  }
  if (any(is.na(rhos))!=TRUE){
  f1<-metacor(rhos,samples);
  print(f1);
  #plot(forest(f1));
  cpval<-metacor(rhos,samples,comb.random=TRUE)$pval.random;
  efs<-metacor(rhos,samples,comb.random=TRUE)$TE.random;
  efsizes<-c(efsizes,efs);
  cpvals<-c(cpvals,cpval);
  }}
  out$efsizes<-efsizes;
  out$cpvals<-cpvals;
  out$V1<-NULL
  return(out)
}

