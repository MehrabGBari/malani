######################################################
######################################################
######################################################
###### Developed by Mehrab Ghanat Bari 3-1-2016#######
## Main code of MALAN method.
## Link to the paper: 
## https://www.nature.com/articles/s41598-017-07481-5

##This code applies  SVM, RandomForest and Generalized linear classifiers
##on Pan-Cancer data (including 9 cancer type) explained on MALAN paper 
##(see link above). The input of this code is the Transcriptomic data from 
##one cancer type and the output is a list of pairs of genes, relevant to
##that cancer. Sice the output of MALANI method is pairs of genes, then 
##the output can be used to construct an inferred network.


##-----------------Required Packages and Directories----------------------
rm(list = ls())
library(class)       # To implement GLM
library(e1071)       # To implement SVM
library(foreach)     # Parralel for loop
library(doSNOW)      # To run code on cluster
library(randomForest)# To implement Random Forest
library(CORElearn)   # To implement 5 different feature selection methods
library(FSelector)   # To implement 5 different feature selection methods
library(xlsx)        # To save results in xlsx format

mainDir <- "/data2/syspharm/projects/03072016_Mehrab_FS"
codeDir <- paste(mainDir,'Codes',sep = '/')
dataDir <- paste(mainDir,'Data',sep = '/')
saveDir <- paste(mainDir,'Results',sep = '/')
source(paste(codeDir,"msvmRFE.R",sep = "/"))

args <- commandArgs()# select "pancreas","prostate","liver",breast" ...
tissues <- args[4]

##-----------------Preprocessing- to convert data in Genes*Samples 
#format with related annotations ----------------------

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

nr <- FALSE#as.logical(args[5])#FALSE#TRUE# normalize each gene expression to N~(0,1)
if(nr){
  expProp <- apply(expProp, 2, function(x){(x-mean(x))/sd(x)}) 
  expProp <- expProp + abs(min(expProp))
}

##-------------------------Parameters------------------------------------
nnodes <- 64                    # Number of nodes for parallel processing
nfolds <- 10                    # k in k-fold cross validation
N = ncol(expProp)               # Number of genes
S <- nrow(expProp)              # Number of samples 
n.Panel <- N-1                  # Number of panels
n.ktop.panels <- floor(0.05*N)  # Top 5% of genes are selected after step 1
gene <- 1:N                     # Gene's index
##----------Divide data to k = nfold for cross validation-----------------
test.id.1 <- rep(1:nfolds,len = S)[sample(S)]
test.id.1 <- lapply(1:nfolds, function(x){which(test.id.1==x)})

##-------------------------Save some deails ------------------------------------
details <- paste(tissues,
                 "N panels length of N-1, top 5 percent panels are selected.",
                 "dot product of of expressions", "normalization", nr, sep = "| "
)
SaveDate <- paste(tissues,"_",format(Sys.time(), "%a%d%m%Y-%H%M"),sep = "")
info.txt <- paste(saveDir,'/info_',SaveDate,'.txt',sep = "")
cat("program has started ...", SaveDate,", dot product, ","normal: ",nr, file = info.txt, sep = "\n")

##------------step 1, N panels to select top 5% of genes -----------------

info <- c()
st.1 <- proc.time()[3]
for(idx.1 in 1:length(tissues)){
  st.2 <- proc.time()[3]
  
  cl <- makeCluster(nnodes,type = 'SOCK')
  registerDoSNOW(cl)
  acc.panels <- foreach(idx.2=1:nfolds,.combine = rbind) %do%{
    
    train.data <- expProp[-test.id.1[[idx.2]],]
    test.data <- expProp[test.id.1[[idx.2]],]
    
    grp.train <- grp.status[-test.id.1[[idx.2]]]
    grp.test <- grp.status[test.id.1[[idx.2]]]
    
    foreach (idx.3=1:N, .combine=cbind, .packages = c("e1071","randomForest")) %dopar% {
      ID <- cbind(rep(idx.3,n.Panel),gene[!gene %in% idx.3])
      tmp.tr <- train.data[,ID[,1]]*train.data[,ID[,2]]
      tmp.ts <- test.data[,ID[,1]]*test.data[,ID[,2]]
      model.svm <- svm(train.data.temp,grp.train,type="C-classification")#,cost=10,scale=F, type="C-classification", kernel="linear")
      svm.pred <- predict(model.svm,test.data.temp)
      ((sum(svm.pred == grp.test)/length(svm.pred))*100)
      
      dat.tr <- data.frame(t(apply(tmp.tr,1,summary)),grp = grp.train)
      dat.ts <- data.frame(t(apply(tmp.ts,1,summary)),grp = grp.test)
      
      m.svm <- svm(dat.tr[,-7],dat.tr[,7],type = "C-classification")
      p.svm <- predict(m.svm,dat.ts[,-7])
      p.svm <- 100*mean(p.svm == dat.ts[,7]) 
      
      m.rf  <- randomForest(grp~., data = dat.tr)
      p.rf <- predict(m.rf,dat.ts,type = "class")
      p.rf <- 100*mean(p.rf == dat.ts$grp)
      
      dat.trLR <- dat.tr;dat.trLR$grp = dat.trLR$grp=="1"
      dat.tsLR <- dat.ts;dat.tsLR$grp = dat.tsLR$grp=="1"
      m.lr <- glm(grp~.,family = binomial,data = dat.trLR)
      p.lr <- predict(m.lr, dat.tsLR[,-7],type = "response")
      p.lr <- 100*sum((p.lr > .5) == dat.tsLR[,7])/nrow(dat.tsLR)
      
      round(c(p.svm,p.rf,p.lr),2)
      
    }
    
  }
  stopCluster(cl)
  
  
  cat(paste("| data:",tissues[idx.1],"| first step done!",
            "| elapsed time  (H) :", (proc.time()[3]-st.1)/3600),
      file = info.txt, append = T,sep = "\n")
  
  dat <- sapply(1:N,function(i) rowMeans(matrix(acc.panels[,i],nrow  = 3,ncol = nfolds)))
  dat <- cbind(1:N,t(dat))
  met <- c("SVM","RF","LR")
  for(idx.2 in 1:length(met)){
    dat.t <- dat[order(dat[,(idx.2+1)],decreasing = T),]
    dat.t <- dat.t[1:n.ktop.panels,c(1,(idx.2+1))];colnames(dat.t) <- c("Panel.id.","Panel.acc.")
    write.xlsx(dat.t, file = paste(saveDir,"/","SVMRFLR_",SaveDate,".xlsx",sep = ""),
               sheetName = met[idx.2],append = T,row.names = F)
  }

  ##------------step 2, combination of 2 of top 5% of genes with all 
  # other genes----------------------------------------------------------
  
  st.1 <- proc.time()[3]
  cl <- makeCluster(nnodes,type = 'SOCK')
  registerDoSNOW(cl)
  met <- c("SVM","RF","LR")
  fn <- paste(saveDir,"/","SVMRFLR_",SaveDate,".xlsx",sep = "")
  for(idx in 1:length(met)){
    top.Panels <- read.xlsx(file = fn, sheetName = met[idx])
    
    top.pairs <- c()
    st.3 <- proc.time()[3]
    for(idx.2 in 1:n.ktop.panels){
      ID <- cbind(rep(top.Panels[idx.2,1],n.Panel),gene[!gene %in% top.Panels[idx.2,1]])
      

      if (met[idx] == "SVM"){
        acc.panels <- foreach(idx.3=1:nfolds,.combine = cbind) %do%{
          train.data <- expProp[-test.id.1[[idx.3]],]
          test.data <- expProp[test.id.1[[idx.3]],]
          
          grp.train <- grp.status[-test.id.1[[idx.3]]]
          grp.test <- grp.status[test.id.1[[idx.3]]]
          
          foreach (idx.4=1:n.Panel, .combine=rbind, .packages = c("e1071","randomForest")) %dopar% {
            
            dat.tr <- data.frame(train.data[,ID[idx.4,]],grp = grp.train)
            dat.ts <- data.frame(test.data[,ID[idx.4,]],grp = grp.test)
            
            m.svm <- svm(dat.tr[,-3],dat.tr[,3],type = "C-classification")
            p.svm <- predict(m.svm,dat.ts[,-3])
            100*mean(p.svm == dat.ts[,3]) 
            
          }
          
        }
      } else if(met[idx] == "RF"){
        acc.panels <- foreach(idx.3=1:nfolds,.combine = cbind) %do%{
          train.data <- expProp[-test.id.1[[idx.3]],]
          test.data <- expProp[test.id.1[[idx.3]],]
          
          grp.train <- grp.status[-test.id.1[[idx.3]]]
          grp.test <- grp.status[test.id.1[[idx.3]]]
          
          foreach (idx.4=1:n.Panel, .combine=rbind, .packages = c("e1071","randomForest")) %dopar% {
            
            dat.tr <- data.frame(train.data[,ID[idx.4,]],grp = grp.train)
            dat.ts <- data.frame(test.data[,ID[idx.4,]],grp = grp.test)
            
            m.rf  <- randomForest(grp~., data = dat.tr)
            p.rf <- predict(m.rf,dat.ts,type = "class")
            100*mean(p.rf == dat.ts$grp)
            
          }
          
        }
        
      }else if(met[idx] == "LR"){
        acc.panels <- foreach(idx.3=1:nfolds,.combine = cbind) %do%{
          train.data <- expProp[-test.id.1[[idx.3]],]
          test.data <- expProp[test.id.1[[idx.3]],]
          
          grp.train <- grp.status[-test.id.1[[idx.3]]]
          grp.test <- grp.status[test.id.1[[idx.3]]]
          
          foreach (idx.4=1:n.Panel, .combine=rbind, .packages = c("e1071","randomForest")) %dopar% {
            
            dat.tr <- data.frame(train.data[,ID[idx.4,]],grp = grp.train)
            dat.ts <- data.frame(test.data[,ID[idx.4,]],grp = grp.test)
            
            dat.trLR <- dat.tr;dat.trLR$grp = dat.trLR$grp=="1"
            dat.tsLR <- dat.ts;dat.tsLR$grp = dat.tsLR$grp=="1"
            m.lr <- glm(grp~.,family = binomial,data = dat.trLR)
            p.lr <- predict(m.lr, dat.tsLR[,-3],type = "response")
            100*sum((p.lr > .5) == dat.tsLR[,3])/nrow(dat.tsLR)
            
          }
          
        }
        
      }
      
      
      dat <- cbind(ID,rowMeans(acc.panels))
      dat <- dat[order(dat[,3],decreasing = T),]
      top.pairs <- rbind(top.pairs,dat[1:10,])

         }
    top.pairs <- cbind(cbind(colnames(expProp)[top.pairs[,1]],colnames(expProp)[top.pairs[,2]]),top.pairs)
    row.names(top.pairs) <- NULL;colnames(top.pairs) <- c("gene.1","gene.2","gene.1.id","gene.2.id","acc.")
    write.xlsx(top.pairs, file = fn,sheetName = paste(met[idx],"step2",sep = "-"),append = T,row.names = F) 
    cat(paste(met[idx]," ","| pair ",idx.2,"out of ",n.ktop.panels, "| elapsed time (H) :", (proc.time()[3]-st.3)/3600,
     "| over all time (H): ",(proc.time()[3]-st.1)/3600), file = info.txt, append = T,sep = "\n")
    
  }
  
  stopCluster(cl)
  
##------------step 3, Apply 5 feature selection methods ------------------
  
  fn <- paste(saveDir,'pancreas_run2.xlsx',sep = "/")
  top.pairs <- read.xlsx(fn,sheetName = "RF-step2")
  
    id <- apply(top.pairs[,c("gene.1.id","gene.2.id")], 2, as.numeric)

  cl <- makeCluster(8,type = "SOCK")
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








