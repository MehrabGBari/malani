######################################################
######################################################
######################################################
###### Developed by Mehrab Ghanat Bari 3-1-2016#######
###### email: m.ghanatbari@gmail.com #################
###### Copyright (C) 2016 ############################
rm(list = ls())
## load required libraries
library(class)
library(e1071)
library(xlsx)
library(foreach)
library(doSNOW)
library(FSelector)
library(CORElearn)
source(paste(codeDir,"msvmRFE.R",sep = "/"))# code for SVM-RFE method developed by johncolby --> http://github.com/johncolby/SVM-RFE
## files and codes location
mainDir <- "/data2/syspharm/projects/03072016_Mehrab_FS"  # main directory
codeDir <- paste(mainDir,'Codes',sep = '/') # code directory
dataDir <- paste(mainDir,'Data',sep = '/') # data directory
saveDir <- paste(mainDir,'Results',sep = '/') # save directory
## set parameters
tissues <- c("breast")# tumor type to be analyzed
#("colon","kidney","liver","lung","ovary","pancreas","prostate","skin)
nnodes <- 64# number of nodes for parallel processing
nr <- FALSE ## map each gene expressions to N~(0,1)

nfolds.1 <- 2
nfolds.2 <- 10
N = ncol(expProp)# number of genes
S <- nrow(expProp)# number of samples
n.Panel <- N-1#
n.ktop.panels <- floor(0.05*N)
gene <- 1:N


if(file.exists(paste(dataDir,"/",tissues,".R",sep = ""))){
  load(paste(dataDir,"/",tissues,".R",sep = ""))
}else{
  load(paste(dataDir,"expProp_PAN_CANCER_hgu133plus2_Mar_21_2016.R",sep = '/'))
  expProp <- t(expProp) #S*N
  load(paste(dataDir,"samplesTable_pan_cancer_microarray_dataset_March_21_2016.R",sep = '/'))

  expProp <- expProp[stAll[,"description1"] == tissues,]
  grp.status <- stAll[stAll[,"description1"] == tissues,"tissue"]
  grp.tissue <- rep(tissues,length(grp.status))
  grp.status <-  c(rep(0,sum(grp.status=="normal")),rep(1,sum(grp.status=="tumor")))
  grp.status <- as.factor(grp.status)
  save(expProp,grp.status,grp.tissue,tissues,file = paste(dataDir,"/",tissues,".R",sep = ""))
}



test.id.1 <- rep(1:nfolds.1,len = S)[sample(S)]
test.id.1 <- lapply(1:nfolds.1, function(x){which(test.id.1==x)})
test.id.2 <- rep(1:nfolds.2,len = S)[sample(S)]
test.id.2 <- lapply(1:nfolds.2, function(x){which(test.id.2==x)})

# if(nr){
#   expProp <- apply(expProp, 2, function(x){(x-mean(x))/sd(x)})
#   expProp <- expProp + abs(min(expProp))
# }

details <- paste(tissues,
                 "N panels length of N-1, top 5 percent panels are selected.",
                 "dot product of expressions","Normalization: ",nr, sep = "| "
)
SaveDate <- paste(tissues,"_",format(Sys.time(), "%a%d%m%Y-%H%M"),sep = "")# paste(tissues,"_",make.names(date()), ".Rout")
info.txt <- paste(saveDir,'/info_',SaveDate,'.txt',sep = "")
cat("program has started ...", SaveDate,", dot product, ","normal: ",nr, file = info.txt, sep = "\n")
#1. 10^6 panels length of 100
####################step 1, N panels ###################################
info <- c()
st.1 <- proc.time()[3]
for(idx.1 in 1:length(tissues)){
  st.2 <- proc.time()[3]

  cl <- makeCluster(nnodes,type = 'SOCK')
  registerDoSNOW(cl)
  acc.pred.panels <- foreach(idx.2=1:nfolds.1,.combine = rbind) %do%{

    train.data <- expProp[-test.id.1[[idx.2]],]
    test.data <- expProp[test.id.1[[idx.2]],]

    grp.train <- grp.status[-test.id.1[[idx.2]]]
    grp.test <- grp.status[test.id.1[[idx.2]]]

    foreach (idx.3=1:N, .combine=cbind, .packages = "e1071") %dopar% {
      ID <- cbind(rep(idx.3,n.Panel),gene[!gene %in% idx.3])
      train.data.temp <- train.data[,ID[,1]]*train.data[,ID[,2]]
      test.data.temp <- test.data[,ID[,1]]*test.data[,ID[,2]]
      model.svm <- svm(train.data.temp,grp.train,type="C-classification")#,cost=10,scale=F, type="C-classification", kernel="linear")
      svm.pred <- predict(model.svm,test.data.temp)
      ((sum(svm.pred == grp.test)/length(svm.pred))*100)
    }

  }
  stopCluster(cl)

  cat(paste("| data:",tissues[idx.1],"| first step done!",
            "| elapsed time  (H) :", (proc.time()[3]-st.1)/3600),
      file = info.txt, append = T,sep = "\n")

  acc.pred.panels <- colMeans(acc.pred.panels)
  acc.pred.panels <- cbind(1:N,acc.pred.panels)
  acc.pred.panels <- acc.pred.panels[order(acc.pred.panels[,2],decreasing = T),]

  top.Panels <- acc.pred.panels[1:n.ktop.panels,]
  top.Panels <- cbind(top.Panels,100-top.Panels[,2])
  colnames(top.Panels) <- c("Panel.id.","Panel.acc.","100_Panel.acc.")

  write.xlsx(top.Panels, file = paste(saveDir,"/","results_",SaveDate,".xlsx",sep = ""),
             sheetName = paste(tissues[idx.1],"_TopPanels",sep = ""),append = T,row.names = F)

  ####################step 2, top 0.05*N panels ###################################
  top.pairs <- c()
  st.3 <- proc.time()[3]
  for(idx.2 in 1:n.ktop.panels){
    ID <- cbind(rep(top.Panels[idx.2,1],n.Panel),gene[!gene %in% top.Panels[idx.2,1]])

    cl <- makeCluster(nnodes,type = 'SOCK')
    registerDoSNOW(cl)
    acc.pred.pairs<- foreach(idx.3=1:nfolds.1,.combine = cbind) %do%{
      train.data <- expProp[-test.id.1[[idx.3]],]
      test.data <- expProp[test.id.1[[idx.3]],]

      grp.train <- grp.status[-test.id.1[[idx.3]]]
      grp.test <- grp.status[test.id.1[[idx.3]]]

      foreach (idx.4=1:n.Panel, .combine=rbind, .packages = "e1071") %dopar% {

        train.data.temp <- train.data[,ID[idx.4,]]
        test.data.temp <- test.data[,ID[idx.4,]]
        model.svm <- svm(train.data.temp,grp.train,type="C-classification")#,cost=10,scale=F, type="C-classification", kernel="linear")
        svm.pred <- predict(model.svm,test.data.temp)
        ((sum(svm.pred == grp.test)/length(svm.pred))*100)
      }

    }
    stopCluster(cl)
    if(idx.2%%50 == 0){
      cat(paste("| pair ",idx.2,"out of ",n.ktop.panels, "| elapsed time (H) :", (proc.time()[3]-st.3)/3600,
                "| over all time (H): ",(proc.time()[3]-st.1)/3600), file = info.txt, append = T,sep = "\n")}

    acc.pred.pairs <- rowMeans(acc.pred.pairs)
    acc.pred.pairs <- cbind(ID,acc.pred.pairs)
    acc.pred.pairs <- cbind(acc.pred.pairs,acc.pred.pairs[,3] - top.Panels[idx.2,2])
    acc.pred.pairs <- acc.pred.pairs[order(acc.pred.pairs[,4],decreasing = T),]
    top.pairs <- rbind(top.pairs,acc.pred.pairs[1:10,])

  }
  cat(paste("| data:",tissues[idx.1],"| Second step done!",
            "| elapsed time  (H) :", (proc.time()[3]-st.3)/3600),
      file = info.txt, append = T,sep = "\n")

  row.names(top.pairs) <- NULL
  top.pairs <- cbind(top.pairs,top.pairs[,1:2])
  colnames(top.pairs) <- c("gene.1","gene.2","acc.","acc.-Panel.acc.","gene.1.id","gene.2.id")
  top.pairs[,1:2] <- cbind(colnames(expProp)[top.pairs[,1]],colnames(expProp)[top.pairs[,2]])

  write.xlsx(top.pairs, file = paste(saveDir,"/","results_",SaveDate,".xlsx",sep = ""),
             sheetName = paste(tissues[idx.1],"PSI",sep = "-"),append = T,row.names = F)

  ####################step 3, FS methods ###################################
  id <- apply(top.pairs[,c("gene.1.id","gene.2.id")], 2, as.numeric)

  cl <- makeCluster(nnodes,type = "SOCK")
  registerDoSNOW(cl)
  feats <- foreach(idx.2=1:nfolds.2,.combine = rbind,.packages = c("FSelector","CORElearn","e1071")) %dopar%{
    train.data <- expProp[-test.id.2[[idx.2]],]
    train.data <- train.data[,id[,1]]+train.data[,id[,2]]

    row.names(train.data) <- NULL
    colnames(train.data) <- NULL

    dimnames(train.data) <- list(rownames(train.data, do.NULL = F, prefix = "sample"),
                                 colnames(train.data, do.NULL = F, prefix = "gene"))
    train.data <- cbind(grp.status[-test.id.2[[idx.2]]],train.data)
    colnames(train.data)[1] <- "grp.status"

    train.data <- as.data.frame(train.data)
    train.data[,-1] <- apply(train.data[,-1], 2, function(x){as.numeric(as.character(x))})

    svm.f <- svmRFE(train.data, k=1, halve.above=100)

    IG.f <- attrEval(grp.status ~ ., train.data, estimator="InfGain")
    IG.f <- sort(IG.f, index.return = T, decreasing = T)$ix

    RF.f <- attrEval(grp.status ~ ., train.data, estimator="ReliefFexpRank", ReliefIterations=30)
    RF.f <- sort(RF.f, index.return = T, decreasing = T)$ix

    SU.f <- symmetrical.uncertainty(grp.status ~ ., train.data)
    SU.f <- sort(SU.f[,1], index.return = T, decreasing = T)$ix

    list(svm.f=svm.f,IG.f=IG.f,RF.f=RF.f,SU.f = SU.f)
  }
  stopCluster(cl)
  nFeat = 1000
  svm.f <- matrix(unlist(feats[,1]),nrow(id),nfolds.2)
  IG.f <- matrix(unlist(feats[,2]),nrow(id),nfolds.2)
  RF.f <- matrix(unlist(feats[,3]),nrow(id),nfolds.2)
  SU.f <- matrix(unlist(feats[,4]),nrow(id),nfolds.2)

  svm.f <- sort(apply(apply(svm.f,2 , function(x){sort(x,index.return=T)$ix}), 1, mean),index.return=T)$ix
  IG.f <- sort(apply(apply(IG.f,2 , function(x){sort(x,index.return=T)$ix}), 1, mean),index.return=T)$ix
  RF.f <- sort(apply(apply(RF.f,2 , function(x){sort(x,index.return=T)$ix}), 1, mean),index.return=T)$ix
  SU.f <- sort(apply(apply(SU.f,2 , function(x){sort(x,index.return=T)$ix}), 1, mean),index.return=T)$ix
  AC.f <- cbind(as.numeric(top.pairs[,"acc."]),1:nrow(top.pairs))
  AC.f <- AC.f[order(AC.f[,1],decreasing = T),2]

  svm.f=svm.f[1:nFeat]
  IG.f=IG.f[1:nFeat]
  RF.f=RF.f[1:nFeat]
  SU.f=SU.f[1:nFeat]
  AC.f=AC.f[1:nFeat]

  comm.f <- cbind(svm.f,IG.f,RF.f,SU.f,AC.f)
  comm.f <- as.data.frame(table(comm.f))
  comm.f <- as.numeric(levels(comm.f[,1]))[comm.f[,1]]

  comm.f <- cbind(comm.f,matrix(0,length(comm.f),6))

  for (idx.3 in 1:nrow(comm.f)){
    f=comm.f[idx.3,1]
    comm.f[idx.3,2] <- sum(c(sum( f == svm.f),sum( f == IG.f),sum( f == RF.f),
                             sum( f == SU.f),sum( f == AC.f)))
    comm.f[idx.3,3:7] <- c(sum( f == svm.f),sum( f == IG.f),sum( f == RF.f),
                           sum( f == SU.f),sum( f == AC.f))
  }
  comm.f <- comm.f[order(comm.f[,2],decreasing = T),]
  comm.f <- cbind(top.pairs[comm.f[,1],c("gene.1","gene.2")],comm.f[,-1])

  colnames(comm.f) <- c("unique.g1","unique.g2","Freq","svm.RFE","InfoGain",
                        "reliefF","SummUn","Acc")
  svm.f <- top.pairs[svm.f,c("gene.1","gene.2")]
  IG.f <- top.pairs[IG.f,c("gene.1","gene.2")]
  RF.f <- top.pairs[RF.f,c("gene.1","gene.2")]
  SU.f <- top.pairs[SU.f,c("gene.1","gene.2")]
  AC.f <- top.pairs[AC.f,c("gene.1","gene.2")]

  results <- cbind(svm.f[1:nFeat,],IG.f[1:nFeat,],RF.f[1:nFeat,],SU.f[1:nFeat,],AC.f[1:nFeat,])
  colnames(results) <- c("svmRFE.g1","svmRFE.g2","InfoGain.g1","InfoGain.g2","reliefF.g1",
                         "reliefF.g2","SummUn.g1","SummUn.g2","Acc.g1","Acc.g2")

  write.xlsx(results, file = paste(saveDir,"/","results_",SaveDate,".xlsx",sep = ""),
             sheetName = "FSmethods",append = T,row.names = F)
  write.xlsx(comm.f, file = paste(saveDir,"/","results_",SaveDate,".xlsx",sep = ""),
             sheetName = "comm.feat",append = T,row.names = F)

  info <- rbind(info,cbind(paste(tissues[idx.1]," run time (H): "),(proc.time()[3]-st.2)/3600))
  info <- rbind(info,cbind("info",details))
  write.xlsx(info, file = paste(saveDir,"/","results_",SaveDate,".xlsx",sep = ""),
             sheetName = "info",append = T, row.names = F,col.names=F)
}








