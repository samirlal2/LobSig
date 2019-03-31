
# Change the following directory
wkdir <- "/Users/samirlal/Desktop/renpjbreastcancernpjbcancer00372majorrevision/"

# Load necessary libraries
library(genefu)
library(caret)
library(matrixStats)
library(ROCR)
library(openxlsx)

setwd(wkdir)
#load("For_Machine_learning.RData")
Combined <- read.table("Combined_for_machine_learning.txt",sep="\t",header=T)
METABRIC <- read.table("METABRIC_for_machine_learning.txt",sep="\t",header=T)
RATHER <- read.table("RATHER_for_machine_learning.txt",sep="\t",header=T)
TCGA <- read.table("TCGA_for_machine_learning.txt",sep="\t",header=T)
Outcome <- read.table("Combined_cohort_classification.txt",sep="\t",header=T)
key.genes<-read.table("14_selected_genes.txt",sep="\t",header=TRUE,row.names=1)
res.sig<-read.table("LobSig.txt",header=TRUE)


Outcome <- Outcome[,c(2,4)]

Outcome$outcome <- ifelse(Outcome$outcome=="not.dead","A","D")

Combined <- merge(Combined,Outcome,by.x="row.names",by.y="Samples")

row.names(Combined) <- Combined$Row.names

Combined$Row.names <- NULL

colnames(Combined)[ncol(Combined)] <- "Class"

####Convert matrix to Z-score



#####Map genes to entrez identifiers
key.genes<-res.sig[res.sig$probe %in% row.names(key.genes),]

####14 gene classifiers
METABRIC_14<-METABRIC[,c(which(colnames(METABRIC) %in% key.genes$probe),ncol(METABRIC))]
RATHER_14<-RATHER[,c(which(colnames(RATHER) %in% key.genes$probe),ncol(RATHER))]
TCGA_14<-TCGA[,c(which(colnames(TCGA) %in% key.genes$probe),ncol(TCGA))]
Combined_14<-Combined[,c(which(colnames(Combined) %in% key.genes$probe),ncol(Combined))] #192 genes


####194 gene classifiers
METABRIC_194<-METABRIC[,c(which(colnames(METABRIC) %in% res.sig$probe),ncol(METABRIC))]
RATHER_194<-RATHER[,c(which(colnames(RATHER) %in% res.sig$probe),ncol(RATHER))]
TCGA_194<-TCGA[,c(which(colnames(TCGA) %in% res.sig$probe),ncol(TCGA))]
Combined_194<-Combined[,c(which(colnames(Combined) %in% res.sig$probe),ncol(Combined))] #192 genes



# each sample rescale the expression values 

datasets <- list(METABRIC_194,RATHER_194,TCGA_194,Combined_194,METABRIC_14,RATHER_14,TCGA_14,Combined_14)


for(i in 1:length(datasets)){
dataset1 <- datasets[[i]]
scale.res<-lapply(1:nrow(dataset1),function(x) {
  y <- dataset1[x,c(1:(ncol(dataset1)-1))] ; 
  new.y <- ((y-min(y))/(max(y)-min(y)))*100
  return(new.y)
})
scaled.res <- as.data.frame(do.call("rbind",scale.res))
scaled.res$Class <- datasets[[i]]$Class

datasets[[i]] <- scaled.res
}

for(i in 1:length(datasets)){
  if(ncol(datasets[[i]])>100){
    datasets[[i]]$LobSig.score <- sig.score(x=res.sig,data=datasets[[i]],annot=res.sig,signed=TRUE)$score 
  }
  else{
    datasets[[i]]$LobSig.score <- sig.score(x=key.genes,data=datasets[[i]],annot=res.sig,signed=TRUE)$score
  }
}


### Re-order the samples 
for(i in 1:length(datasets)){
  datasets[[i]] <- datasets[[i]][,c(1:(ncol(datasets[[i]])-2),ncol(datasets[[i]]),(ncol(datasets[[i]])-1))]
  }



###single gene model 


training_results<-list()
for(n in 1:length(datasets)){
# Apply 5 k-cv in R 
datas<-datasets[[n]]
#Create 5 equally size folds
folds<-createFolds(datas$Class, k = 5, list = TRUE, returnTrain = FALSE)
all_folds<- list()
#Perform 5 fold cross validation
for(i in 1:5){
  #Segement your data by fold function 
  testIndexes <- folds[[i]]
  testData <- datas[testIndexes, ]
  trainData <- datas[-testIndexes, ]
  Class <- as.vector(trainData[,ncol(trainData)])
  Score <- as.vector(as.numeric(trainData[,(ncol(trainData)-1)]))
  #Use the test and train data partitions however you desire...
  pred<-ROCR::prediction(Score,Class)
  ss <- performance(pred, "sens", "spec")
  cutoff <- ss@alpha.values[[1]][which.max(ss@x.values[[1]]+ss@y.values[[1]])]
  #Perform classifications here 
  testData[testData$LobSig.score<cutoff,(ncol(testData)+1)]<-"Low"
  testData[testData$LobSig.score>=cutoff,ncol(testData)]<-"High"
  all_folds[[i]] <- testData
}
all_folds <- do.call("rbind", all_folds)
training_results[[n]] <- all_folds
}

training_results <- training_results[1:4]

names(training_results) <- c("METABRIC_194","RATHER_194","TCGA_194","Combined_194")

write.xlsx(training_results, file = "Combined_cohort_classification.xlsx")
