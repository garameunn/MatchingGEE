# library(Matrix, lib.loc="/home2/nekim/Rpackage/")
rm(list=ls())
options(digits = 3)
# .libPaths("/data/pkg36")
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/DealFunc.R")
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") # 
# load(file="/data4/nekim/matching/Rraw/qcotu_toy.RData") # otutable1 = filtered qc otu
# seeds = sample.int(10E+6, replace=FALSE);
# save(seeds, file="/home2/nekim/scratch/matching/Robject/seed.RData")
load(file="/home2/nekim/scratch/matching/Robject_old/seed.RData")
bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

load(file="/home2/nekim/scratch/matching/Robject/otulist.RData") # otulist (otu (1000/500/200), lbsize (1000))
load(file="/home2/nekim/scratch/matching/Robject/indiclist.RData") # indiclist
load(file="/home2/nekim/scratch/matching/Robject/psdifflist.RData") # psdifflist
nsamchar <- c(1000,500,200) # norder=1





#####################################
## 챕터 : 2 (Binary)
## 시나리오: 1
## 통계량 : T1
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : SourceCode_Binary_Scene123.R
#####################################

#### nsamp 변경
nsamp = 200

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = 0

scenario="S1"

#### rate, rep, nonlin 변경
set.seed(1)
rep=5000
dimm=dim(otutable)[1]; omics="Metagenomics"
datatable=otutable;rate=0.5;nonlin="comp";
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    tablefun <- get("tablemaking_t1")
    tablefun1 <- get("tablemaking_t1_wo")
    
    
    tble0 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.05 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.1 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.2 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.3 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=NULL, nonlin=nonlin)
    
    tble0_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=NULL, nonlin=nonlin) 
    tble0_.05_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.1_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.2_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.3_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=NULL, nonlin=nonlin)
    
    
    # 
    r.sig <- coeff <- r.sig2 <- coeff2 <- c() # withrep nocal (Indep - EX)
    r.sig0.1 <- coeff0.1 <- r.sig20.1 <- coeff20.1 <- c() # withrep cal0.1 
    r.sig0.2 <- coeff0.2 <- r.sig20.2 <- coeff20.2 <- c() # withrep cal0.2 
    r.sig0.3 <- coeff0.3 <- r.sig20.3 <- coeff20.3 <- c() # withrep cal0.3
    r.sig0.05 <- coeff0.05 <- r.sig20.05 <- coeff20.05 <- c() # withrep cal0.05
    
    r.sig_w <- coeff_w <- r.sig2_w <- coeff2_w <- c() # without rep nocal (Indep - EX)
    r.sig0.1_w <- coeff0.1_w <- r.sig20.1_w <- coeff20.1_w <- c() # withrep cal0.1 
    r.sig0.2_w <- coeff0.2_w <- r.sig20.2_w <- coeff20.2_w <- c() # withrep cal0.2 
    r.sig0.3_w <- coeff0.3_w <- r.sig20.3_w <- coeff20.3_w <- c() # withrep cal0.3
    r.sig0.05_w <- coeff0.05_w <- r.sig20.05_w <- coeff20.05_w <- c() # withrep cal0.05
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- c()
    
  
    tble <- data.frame(tble0[[1]])
    tble0.1 <- tble0_.1[[1]]
    tble0.2 <- tble0_.2[[1]]
    tble0.3 <- tble0_.3[[1]]
    tble0.05 <- tble0_.05[[1]]
    tble_wo <- tble0_wo[[1]]
    tble0.1_wo <- tble0_.1_wo[[1]]
    tble0.2_wo <- tble0_.2_wo[[1]]
    tble0.3_wo <- tble0_.3_wo[[1]]
    tble0.05_wo <- tble0_.05_wo[[1]]
    
    tble2 <- data.frame(tble0[[2]]) # entire data
    tble2$cmbt <- unlist(cmbatable[x,])
    
    
    #######################
    ###### with rep #######
    #######################
    # no caliper
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig <- c(r.sig,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff <- c(coeff,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig2 <- c(r.sig2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff2 <- c(coeff2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2); 
    
    # 0.2 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.2 <- c(r.sig0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.2 <- c(coeff0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.2 <- c(r.sig20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.2 <- c(coeff20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    rm(geemfit);rm(geemfit2); 
    
    # 0.1 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.1 <- c(r.sig0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.1 <- c(coeff0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.1 <- c(r.sig20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.1 <- c(coeff20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    rm(geemfit);rm(geemfit2);
    
    # 0.3 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.3 <- c(r.sig0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.3 <- c(coeff0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.3 <- c(r.sig20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.3 <- c(coeff20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.05 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.05 <- c(r.sig0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.05 <- c(coeff0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.05 <- c(r.sig20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.05 <- c(coeff20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    #########################
    ###### without rep ######
    #########################
    rm(geemfit);rm(geemfit2); 
    # no caliper
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig_w <- c(r.sig_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff_w <- c(coeff_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig2_w <- c(r.sig2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff2_w <- c(coeff2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.2 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.2_w <- c(r.sig0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.2_w <- c(coeff0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.2_w <- c(r.sig20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.2_w <- c(coeff20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.1 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.1_w <- c(r.sig0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.1_w <- c(coeff0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.1_w <- c(r.sig20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.1_w <- c(coeff20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    # 0.3 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.3_w <- c(r.sig0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.3_w <- c(coeff0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.3_w <- c(r.sig20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.3_w <- c(coeff20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.05 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.05_w <- c(r.sig0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.05_w <- c(coeff0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.05_w <- c(r.sig20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.05_w <- c(coeff20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2);
    
    
    logit <- glm(paste("case~x1+PS","marker",sep="+"), data=tble2, family = "binomial")
    logit.sig <- c(logit.sig,ifelse(logit$converged,coef(summary(logit))[4,4],1))
    logit.coef <- c(logit.coef,ifelse(logit$converged,coef(summary(logit))[4,1],NA))
    rm(logit)
    
    logit1 <- glm(paste("case~x1","marker",sep="+"), data=tble2, family = "binomial")
    logit.sig1 <- ifelse(logit1$converged,coef(summary(logit1))[3,4],1)
    logit.coef1 <- ifelse(logit1$converged,coef(summary(logit1))[3,1],NA)
    rm(logit1)
    
    
    combat <- glm(paste("case~x1+hetero1+hetero2","cmbt",sep="+"), data=tble2, family = "binomial") 
    combat.sig <- c(combat.sig,ifelse(combat$converged,coef(summary(combat))[5,4],1))
    combat.coef <- c(combat.coef,ifelse(combat$converged,coef(summary(combat))[5,1],NA))
    
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, 
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #23
                      logit.coef1, logit.coef, combat.coef, 
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w, #23
                      nrow(tble), nrow(tble0.3), nrow(tble0.2), nrow(tble0.1), nrow(tble0.05),
                      nrow(tble_wo), nrow(tble0.3_wo), nrow(tble0.2_wo), nrow(tble0.1_wo), nrow(tble0.05_wo), #10
                      tble0[[3]], tble0_.3[[3]], tble0_.2[[3]], tble0_.1[[3]], tble0_.05[[3]])) #5 (meanrep)
    
    
  },mc.cores=6))
},mc.cores = 6)

pvalues_comb <- do.call(rbind,lapply(pvalues,function(x) matrix(unlist(x), byrow=T,ncol=61))) # 두번째 mclapply 앞에 unlist 씌우는 대신 power 와 동일하게 씌우지 않고, 정리할 때 씌움
head(pvalues_comb)

save(pvalues_comb, file=paste0("/data4/nekim/matching/Robject/T1/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 
# load("/data4/nekim/matching/Robject/T1/NSAMP1000_REP5000_RATE0.2_inv_S1.RData")

t1err <- apply(pvalues_comb[,1:23], 2, type1_1)
rownames(t1err) <-c("p.001","p.005","p.01","p.05","p.1")
t1err
CIfun1<-function(pvals, level){ # 틀린 CI
  # browser()
  if(sum(is.na(pvals))>0) pvals[which(is.na(pvals))] <- 1
  xx <-round((sum(pvals<level)/length(pvals))*5000,0)
  testp <- prop.test(x=xx,n=5000,correct=T)
  return(testp$conf.int)
}
conf.int <- do.call(rbind,lapply(c(0.001,0.005,0.01,0.05,0.1),function(y) apply(pvalues_comb[,1:22], 2, FUN=CIfun1, level=y))) # conf.int
rownames(conf.int) <- paste0("CI_",c(rep(c("p.001","p.005","p.01","p.05","p.1"),each=2)))
t1err1 <- rbind(t1err, conf.int)
colnames(t1err1) <-  c("Unadjusted",wghtMnames)
t1err1
t1err2 <- t1err1
t(t1err2)





#####################################
## 챕터 : 2 (Binary)
## 시나리오: 2
## 통계량 : T1
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : 
#####################################

#### nsamp 변경
nsamp = 200

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = 0

scenario="S2"


#### rate, rep, nonlin 변경
set.seed(1)
rep=5000
dimm=dim(otutable)[1]; omics="Metagenomics"
datatable=otutable;rate=0.2;nonlin=NULL; # S1 이 아니면 nonlin 은 NULL
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    tablefun <- get("tablemaking_t1")
    tablefun1 <- get("tablemaking_t1_wo")
    
    
    tble0 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.05 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.1 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.2 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.3 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=NULL, nonlin=nonlin)
    
    tble0_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=NULL, nonlin=nonlin) 
    tble0_.05_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.1_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.2_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.3_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=NULL, nonlin=nonlin)
    
    
    # 
    r.sig <- coeff <- r.sig2 <- coeff2 <- c() # withrep nocal (Indep - EX)
    r.sig0.1 <- coeff0.1 <- r.sig20.1 <- coeff20.1 <- c() # withrep cal0.1 
    r.sig0.2 <- coeff0.2 <- r.sig20.2 <- coeff20.2 <- c() # withrep cal0.2 
    r.sig0.3 <- coeff0.3 <- r.sig20.3 <- coeff20.3 <- c() # withrep cal0.3
    r.sig0.05 <- coeff0.05 <- r.sig20.05 <- coeff20.05 <- c() # withrep cal0.05
    
    r.sig_w <- coeff_w <- r.sig2_w <- coeff2_w <- c() # without rep nocal (Indep - EX)
    r.sig0.1_w <- coeff0.1_w <- r.sig20.1_w <- coeff20.1_w <- c() # withrep cal0.1 
    r.sig0.2_w <- coeff0.2_w <- r.sig20.2_w <- coeff20.2_w <- c() # withrep cal0.2 
    r.sig0.3_w <- coeff0.3_w <- r.sig20.3_w <- coeff20.3_w <- c() # withrep cal0.3
    r.sig0.05_w <- coeff0.05_w <- r.sig20.05_w <- coeff20.05_w <- c() # withrep cal0.05
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- c()
    
    
    tble <- data.frame(tble0[[1]])
    tble0.1 <- tble0_.1[[1]]
    tble0.2 <- tble0_.2[[1]]
    tble0.3 <- tble0_.3[[1]]
    tble0.05 <- tble0_.05[[1]]
    tble_wo <- tble0_wo[[1]]
    tble0.1_wo <- tble0_.1_wo[[1]]
    tble0.2_wo <- tble0_.2_wo[[1]]
    tble0.3_wo <- tble0_.3_wo[[1]]
    tble0.05_wo <- tble0_.05_wo[[1]]
    
    tble2 <- data.frame(tble0[[2]]) # entire data
    tble2$cmbt <- unlist(cmbatable[x,])
    
    
    #######################
    ###### with rep #######
    #######################
    # no caliper
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig <- c(r.sig,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff <- c(coeff,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig2 <- c(r.sig2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff2 <- c(coeff2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2); 
    
    # 0.2 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.2 <- c(r.sig0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.2 <- c(coeff0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.2 <- c(r.sig20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.2 <- c(coeff20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    rm(geemfit);rm(geemfit2); 
    
    # 0.1 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.1 <- c(r.sig0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.1 <- c(coeff0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.1 <- c(r.sig20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.1 <- c(coeff20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    rm(geemfit);rm(geemfit2);
    
    # 0.3 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.3 <- c(r.sig0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.3 <- c(coeff0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.3 <- c(r.sig20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.3 <- c(coeff20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.05 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.05 <- c(r.sig0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.05 <- c(coeff0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.05 <- c(r.sig20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.05 <- c(coeff20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    #########################
    ###### without rep ######
    #########################
    rm(geemfit);rm(geemfit2); 
    # no caliper
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig_w <- c(r.sig_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff_w <- c(coeff_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig2_w <- c(r.sig2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff2_w <- c(coeff2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.2 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.2_w <- c(r.sig0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.2_w <- c(coeff0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.2_w <- c(r.sig20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.2_w <- c(coeff20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.1 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.1_w <- c(r.sig0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.1_w <- c(coeff0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.1_w <- c(r.sig20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.1_w <- c(coeff20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    # 0.3 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.3_w <- c(r.sig0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.3_w <- c(coeff0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.3_w <- c(r.sig20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.3_w <- c(coeff20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.05 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.05_w <- c(r.sig0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.05_w <- c(coeff0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.05_w <- c(r.sig20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.05_w <- c(coeff20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2);
    
    
    logit <- glm(paste("case~x1+PS","marker",sep="+"), data=tble2, family = "binomial")
    logit.sig <- c(logit.sig,ifelse(logit$converged,coef(summary(logit))[4,4],1))
    logit.coef <- c(logit.coef,ifelse(logit$converged,coef(summary(logit))[4,1],NA))
    rm(logit)
    
    logit1 <- glm(paste("case~x1","marker",sep="+"), data=tble2, family = "binomial")
    logit.sig1 <- ifelse(logit1$converged,coef(summary(logit1))[3,4],1)
    logit.coef1 <- ifelse(logit1$converged,coef(summary(logit1))[3,1],NA)
    rm(logit1)
    
    
    combat <- glm(paste("case~x1+hetero1+hetero2","cmbt",sep="+"), data=tble2, family = "binomial") 
    combat.sig <- c(combat.sig,ifelse(combat$converged,coef(summary(combat))[5,4],1))
    combat.coef <- c(combat.coef,ifelse(combat$converged,coef(summary(combat))[5,1],NA))
    
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, 
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #23
                      logit.coef1, logit.coef, combat.coef, 
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w, #23
                      nrow(tble), nrow(tble0.3), nrow(tble0.2), nrow(tble0.1), nrow(tble0.05),
                      nrow(tble_wo), nrow(tble0.3_wo), nrow(tble0.2_wo), nrow(tble0.1_wo), nrow(tble0.05_wo), #10
                      tble0[[3]], tble0_.3[[3]], tble0_.2[[3]], tble0_.1[[3]], tble0_.05[[3]])) #5 (meanrep)
    
    
  },mc.cores=6))
},mc.cores = 6)

pvalues_comb <- do.call(rbind,lapply(pvalues,function(x) matrix(unlist(x), byrow=T,ncol=61))) # 두번째 mclapply 앞에 unlist 씌우는 대신 power 와 동일하게 씌우지 않고, 정리할 때 씌움
head(pvalues_comb)

save(pvalues_comb, file=paste0("/data4/nekim/matching/Robject/T1/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 




#####################################
## 챕터 : 2 (Binary)
## 시나리오: 3
## 통계량 : T1
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : OmicGee_ver14.R
#####################################

#### nsamp 변경
nsamp = 500

targetnum <- which(nsamchar==nsamp)

# otutable <- otulist[[targetnum]]
otutable <- as(otu_table(otulist[[targetnum]]),"matrix") # error in validObject(.Object) : invalid class “otu_table” object:Non-numeric matrix provided as OTU table. 해결
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = 0

scenario="S3"


# k 설정
k <- c(1,8,61,77)



#### rate, rep, nonlin 변경
set.seed(1)
rep=5000
dimm=dim(otutable)[1]; omics="Metagenomics"
datatable=otutable;rate=0.2;nonlin=NULL; # S1 이 아니면 nonlin 은 NULL
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    tablefun <- get("tablemaking_t1")
    tablefun1 <- get("tablemaking_t1_wo")
    
    
    tble0 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=k, nonlin=nonlin)  
    tble0_.05 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=k, nonlin=nonlin)
    tble0_.1 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=k, nonlin=nonlin)
    tble0_.2 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=k, nonlin=nonlin)
    tble0_.3 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=k, nonlin=nonlin)
    
    tble0_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=k, nonlin=nonlin) 
    tble0_.05_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=k, nonlin=nonlin)  
    tble0_.1_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=k, nonlin=nonlin)
    tble0_.2_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=k, nonlin=nonlin)  
    tble0_.3_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=k, nonlin=nonlin)
    
    
    # 
    r.sig <- coeff <- r.sig2 <- coeff2 <- c() # withrep nocal (Indep - EX)
    r.sig0.1 <- coeff0.1 <- r.sig20.1 <- coeff20.1 <- c() # withrep cal0.1 
    r.sig0.2 <- coeff0.2 <- r.sig20.2 <- coeff20.2 <- c() # withrep cal0.2 
    r.sig0.3 <- coeff0.3 <- r.sig20.3 <- coeff20.3 <- c() # withrep cal0.3
    r.sig0.05 <- coeff0.05 <- r.sig20.05 <- coeff20.05 <- c() # withrep cal0.05
    
    r.sig_w <- coeff_w <- r.sig2_w <- coeff2_w <- c() # without rep nocal (Indep - EX)
    r.sig0.1_w <- coeff0.1_w <- r.sig20.1_w <- coeff20.1_w <- c() # withrep cal0.1 
    r.sig0.2_w <- coeff0.2_w <- r.sig20.2_w <- coeff20.2_w <- c() # withrep cal0.2 
    r.sig0.3_w <- coeff0.3_w <- r.sig20.3_w <- coeff20.3_w <- c() # withrep cal0.3
    r.sig0.05_w <- coeff0.05_w <- r.sig20.05_w <- coeff20.05_w <- c() # withrep cal0.05
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- c()
    
    
    tble <- data.frame(tble0[[1]])
    tble0.1 <- tble0_.1[[1]]
    tble0.2 <- tble0_.2[[1]]
    tble0.3 <- tble0_.3[[1]]
    tble0.05 <- tble0_.05[[1]]
    tble_wo <- tble0_wo[[1]]
    tble0.1_wo <- tble0_.1_wo[[1]]
    tble0.2_wo <- tble0_.2_wo[[1]]
    tble0.3_wo <- tble0_.3_wo[[1]]
    tble0.05_wo <- tble0_.05_wo[[1]]
    
    tble2 <- data.frame(tble0[[2]]) # entire data
    tble2$cmbt <- unlist(cmbatable[x,])
    
    
    #######################
    ###### with rep #######
    #######################
    # no caliper
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig <- c(r.sig,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff <- c(coeff,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig2 <- c(r.sig2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff2 <- c(coeff2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2); 
    
    # 0.2 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.2 <- c(r.sig0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.2 <- c(coeff0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.2 <- c(r.sig20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.2 <- c(coeff20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    rm(geemfit);rm(geemfit2); 
    
    # 0.1 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.1 <- c(r.sig0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.1 <- c(coeff0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.1 <- c(r.sig20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.1 <- c(coeff20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    rm(geemfit);rm(geemfit2);
    
    # 0.3 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.3 <- c(r.sig0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.3 <- c(coeff0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.3 <- c(r.sig20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.3 <- c(coeff20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.05 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.05 <- c(r.sig0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.05 <- c(coeff0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.05 <- c(r.sig20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.05 <- c(coeff20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    #########################
    ###### without rep ######
    #########################
    rm(geemfit);rm(geemfit2); 
    # no caliper
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig_w <- c(r.sig_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff_w <- c(coeff_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig2_w <- c(r.sig2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff2_w <- c(coeff2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.2 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.2_w <- c(r.sig0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.2_w <- c(coeff0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.2_w <- c(r.sig20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.2_w <- c(coeff20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.1 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.1_w <- c(r.sig0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.1_w <- c(coeff0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.1_w <- c(r.sig20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.1_w <- c(coeff20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    # 0.3 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.3_w <- c(r.sig0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.3_w <- c(coeff0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.3_w <- c(r.sig20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.3_w <- c(coeff20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    
    
    rm(geemfit);rm(geemfit2);
    
    # 0.05 caliper 
    geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                  corstr = "independence", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig0.05_w <- c(r.sig0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
    coeff0.05_w <- c(coeff0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
    
    geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                   corstr = "exchangeable", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
    
    r.sig20.05_w <- c(r.sig20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
    coeff20.05_w <- c(coeff20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
    
    
    rm(geemfit);rm(geemfit2);
    
    
    logit <- glm(paste("case~x1+PS","marker",sep="+"), data=tble2, family = "binomial")
    logit.sig <- c(logit.sig,ifelse(logit$converged,coef(summary(logit))[4,4],1))
    logit.coef <- c(logit.coef,ifelse(logit$converged,coef(summary(logit))[4,1],NA))
    rm(logit)
    
    logit1 <- glm(paste("case~x1","marker",sep="+"), data=tble2, family = "binomial")
    logit.sig1 <- ifelse(logit1$converged,coef(summary(logit1))[3,4],1)
    logit.coef1 <- ifelse(logit1$converged,coef(summary(logit1))[3,1],NA)
    rm(logit1)
    
    
    combat <- glm(paste("case~x1+hetero1+hetero2","cmbt",sep="+"), data=tble2, family = "binomial") 
    combat.sig <- c(combat.sig,ifelse(combat$converged,coef(summary(combat))[5,4],1))
    combat.coef <- c(combat.coef,ifelse(combat$converged,coef(summary(combat))[5,1],NA))
    
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, 
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #23
                      logit.coef1, logit.coef, combat.coef, 
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w, #23
                      nrow(tble), nrow(tble0.3), nrow(tble0.2), nrow(tble0.1), nrow(tble0.05),
                      nrow(tble_wo), nrow(tble0.3_wo), nrow(tble0.2_wo), nrow(tble0.1_wo), nrow(tble0.05_wo), #10
                      tble0[[3]], tble0_.3[[3]], tble0_.2[[3]], tble0_.1[[3]], tble0_.05[[3]])) #5 (meanrep)
    
    
  },mc.cores=6))
},mc.cores = 10)

pvalues_comb <- do.call(rbind,lapply(pvalues,function(x) matrix(unlist(x), byrow=T,ncol=61))) # 두번째 mclapply 앞에 unlist 씌우는 대신 power 와 동일하게 씌우지 않고, 정리할 때 씌움
head(pvalues_comb)

save(pvalues_comb, file=paste0("/data4/nekim/matching/Robject/T1/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 


















#####################################
#####################################
## 챕터 : 1,2 
## 시나리오: 1,2,3
## 통계량 : T1
## 주제 : 정리 (계산/정리 중 하나)
## 출처 : tmp_T1.R
#####################################
source("/data4/nekim/matching/Rcode/Rcode_DegreeTotal/DealFunc.R")
load(file="/data4/nekim/matching/Robject/EffectSizeList.RData") # b2list 


# nsamp<-500;scenario<-"S1";Rate="0.3";rep=5000;nonlin="inv"
# 
# 
# nsample <- c(200,500,1000)
# scen <- c("S1","S2") # c("S1","S2","S3)
# rate <- c("0.2","0.3","0.5")
# # nonlinear <- ifelse(scenario=="S1",c("inv","comp"),NULL)
# rep=5000
# 
# # 1~23: Pvalue, 24~46: coefficient, 47~56: nrow, 57~61: mean rep
# load(file=paste0("/data4/nekim/matching/Robject/T1/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 
# 
# t1err <- apply(pvalues_comb[,1:23], 2, type1_1)
# rownames(t1err) <-c("p.001","p.005","p.01","p.05","p.1")
# t1err
# CIfun1<-function(pvals, level){ # 틀린 CI
#   # browser()
#   if(sum(is.na(pvals))>0) pvals[which(is.na(pvals))] <- 1
#   xx <-round((sum(pvals<level)/length(pvals))*5000,0)
#   testp <- prop.test(x=xx,n=5000,correct=T)
#   return(testp$conf.int)
# }
# conf.int <- do.call(rbind,lapply(c(0.001,0.005,0.01,0.05,0.1),function(y) apply(pvalues_comb[,1:23], 2, FUN=CIfun1, level=y))) # conf.int
# rownames(conf.int) <- paste0("CI_",c(rep(c("p.001","p.005","p.01","p.05","p.1"),each=2)))
# t1err1 <- rbind(t1err, conf.int)
# colnames(t1err1) <-  c("LogitZ","LogitZX","Combat",
#                        paste0("Cal", rep(c("No",0.3,0.2,0.1,0.05), each=2), c("Ind","Exp")),
#                        paste0("Cal", rep(c("No",0.3,0.2,0.1,0.05), each=2), c("Ind","Exp") ,"Wrep"))
# t1err1






CIfun1<-function(pvals, level){ # 틀린 CI
  # browser()
  if(sum(is.na(pvals))>0) pvals[which(is.na(pvals))] <- 1
  xx <-round((sum(pvals<level)/length(pvals))*5000,0)
  testp <- prop.test(x=xx,n=5000,correct=T)
  return(testp$conf.int)
}

nsample <- c(200, 500, 1000)
scen <- c("S1", "S2") # 나중에 "S3" 추가 가능
rate <- c("0.2", "0.3", "0.5")
rep <- 5000

result_list <- list()

for (nsamp in nsample) {
  for (scenario in scen) {
    if (scenario == "S1") {
      nonlin_values <- c("inv", "comp")
    } else {
      nonlin_values <- NULL
    }
    
    for (rate_value in rate) {
      for (nonlin in if (is.null(nonlin_values)) "" else nonlin_values) {
        file_name <- paste0("NSAMP", nsamp, "_REP", rep, "_RATE", rate_value, "_",
                            if (nonlin != "") nonlin, "_", scenario, ".RData")
        print(file_name)
        load(file = paste0("/data4/nekim/matching/Robject/T1/", file_name))
        
        t1err <- apply(pvalues_comb[,1:23], 2, type1_1)
        rownames(t1err) <- c("p.001","p.005","p.01","p.05","p.1")
        
        conf.int <- do.call(rbind, lapply(c(0.001, 0.005, 0.01, 0.05, 0.1), 
                                          function(y) apply(pvalues_comb[,1:23], 2, FUN = CIfun1, level = y)))
        rownames(conf.int) <- paste0("CI_", c(rep(c("p.001", "p.005", "p.01", "p.05", "p.1"), each = 2)))
        
        t1err1 <- rbind(t1err, conf.int)
        colnames(t1err1) <- c("LogitZ", "LogitZX", "Combat",
                              paste0("Cal", rep(c("No", 0.3, 0.2, 0.1, 0.05), each = 2), c("Ind", "Exp")),
                              paste0("Cal", rep(c("No", 0.3, 0.2, 0.1, 0.05), each = 2), c("Ind", "Exp"), "Wrep"))
        
        result_list[[file_name]] <- t1err1
      }
    }
  }
}

result_list
combined_table <- do.call(rbind, lapply(names(result_list), function(table_name) {
  data <- result_list[[table_name]]
  data <- as.data.frame(data[-c(1,6,7),]) # 데이터프레임으로 변환
  table_name1 <- gsub(pattern = ".RData","",table_name)
  data1 <- data.frame(table_name1, data) # 새로운 열 추가
  return(data1)
}))

head(combined_table)
combined_table <- data.frame(title=rep(rownames(result_list[[1]])[-c(1,6,7)],length(result_list)),
                             combined_table)

write.xlsx(combined_table, file = "/data4/nekim/matching/Routput_240818/combined_T1err_woS3.xlsx", row.names = FALSE)


#### S3 도 포함하는 코드

CIfun1<-function(pvals, level){ # 틀린 CI
  # browser()
  if(sum(is.na(pvals))>0) pvals[which(is.na(pvals))] <- 1
  xx <-round((sum(pvals<level)/length(pvals))*5000,0)
  testp <- prop.test(x=xx,n=5000,correct=T)
  return(testp$conf.int)
}

nsample <- c(200, 500, 1000)
scen <- c("S1", "S2","S3") # 나중에 "S3" 추가 가능
rate <- c("0.2", "0.3", "0.5")
rep <- 5000

result_list <- list()

for (nsamp in nsample) {
  for (scenario in scen) {
    if (scenario == "S1") {
      nonlin_values <- c("inv", "comp")
    } else {
      nonlin_values <- NULL
    }
    
    for (rate_value in rate) {
      for (nonlin in if (is.null(nonlin_values)) "" else nonlin_values) {
        file_name <- paste0("NSAMP", nsamp, "_REP", rep, "_RATE", rate_value, "_",
                            if (nonlin != "") nonlin, "_", scenario, ".RData")
        print(file_name)
        load(file = paste0("/data4/nekim/matching/Robject/T1/", file_name))
        
        t1err <- apply(pvalues_comb[,1:23], 2, type1_1)
        rownames(t1err) <- c("p.001","p.005","p.01","p.05","p.1")
        
        conf.int <- do.call(rbind, lapply(c(0.001, 0.005, 0.01, 0.05, 0.1), 
                                          function(y) apply(pvalues_comb[,1:23], 2, FUN = CIfun1, level = y)))
        rownames(conf.int) <- paste0("CI_", c(rep(c("p.001", "p.005", "p.01", "p.05", "p.1"), each = 2)))
        
        t1err1 <- rbind(t1err, conf.int)
        colnames(t1err1) <- c("LogitZ", "LogitZX", "Combat",
                              paste0("Cal", rep(c("No", 0.3, 0.2, 0.1, 0.05), each = 2), c("Ind", "Exp")),
                              paste0("Cal", rep(c("No", 0.3, 0.2, 0.1, 0.05), each = 2), c("Ind", "Exp"), "Wrep"))
        
        result_list[[file_name]] <- t1err1
      }
    }
  }
}

result_list
combined_table <- do.call(rbind, lapply(names(result_list), function(table_name) {
  data <- result_list[[table_name]]
  data <- as.data.frame(data[-c(1,6,7),]) # 데이터프레임으로 변환
  table_name1 <- gsub(pattern = ".RData","",table_name)
  data1 <- data.frame(table_name1, data) # 새로운 열 추가
  return(data1)
}))

head(combined_table)
combined_table <- data.frame(title=rep(rownames(result_list[[1]])[-c(1,6,7)],length(result_list)),
                             combined_table)

write.xlsx(combined_table, file = "/data4/nekim/matching/Routput_240818/combined_T1err.xlsx", row.names = FALSE)




#####################################
#####################################
## 챕터 : 1 
## 시나리오: 1,2,3
## 통계량 : T1_nrow, mean rep
## 주제 : 정리 (계산/정리 중 하나)
## 출처 : tmp_T1.R
#####################################

nsample <- c(200, 500, 1000)
scen <- c("S1", "S2","S3") # 나중에 "S3" 추가 가능
rate <- c("0.2", "0.3", "0.5")
rep <- 5000

mean_nsample_list <- list()
mean_rep_list <- list()

for (nsamp in nsample) {
  for (scenario in scen) {
    if (scenario == "S1") {
      nonlin_values <- c("inv", "comp")
    } else {
      nonlin_values <- NULL
    }
    
    for (rate_value in rate) {
      for (nonlin in if (is.null(nonlin_values)) "" else nonlin_values) {
        file_name <- paste0("NSAMP", nsamp, "_REP", rep, "_RATE", rate_value, "_",
                            if (nonlin != "") nonlin, "_", scenario, ".RData")
        print(file_name)
        load(file = paste0("/data4/nekim/matching/Robject/T1/", file_name))
        
        # 47~56번째 열 평균값 계산 (CalNo, Cal0.3, Cal0.2, Cal0.1, Cal0.05, CalNoRep, Cal0.3Rep, Cal0.2Rep, Cal0.1Rep, Cal0.05Rep)
        mean_nsample_values <- colMeans(pvalues_comb[, 47:56], na.rm = TRUE)
        names(mean_nsample_values) <- c("CalNo", "Cal0.3", "Cal0.2", "Cal0.1", "Cal0.05", 
                                        "CalNoRep", "Cal0.3Rep", "Cal0.2Rep", "Cal0.1Rep", "Cal0.05Rep")
        
        # mean_nsample 테이블에 파일명과 함께 추가
        mean_nsample_list[[file_name]] <- mean_nsample_values
        
        # 57~61번째 열 평균값 계산 (CalNoRep, Cal0.3Rep, Cal0.2Rep, Cal0.1Rep, Cal0.05Rep)
        mean_rep_values <- colMeans(pvalues_comb[, 57:61], na.rm = TRUE)
        names(mean_rep_values) <- c("CalNoRep", "Cal0.3Rep", "Cal0.2Rep", "Cal0.1Rep", "Cal0.05Rep")
        
        # mean_rep 테이블에 파일명과 함께 추가
        mean_rep_list[[file_name]] <- mean_rep_values
      }
    }
  }
}

# mean_nsample 테이블 결합 (Wide form으로 열로 저장)
mean_nsample <- do.call(rbind, lapply(names(mean_nsample_list), function(table_name) {
  data <- mean_nsample_list[[table_name]]
  table_name1 <- gsub(pattern = ".RData", "", table_name)
  data1 <- as.data.frame(matrix(data, ncol = 10, byrow = T)) # 값을 행으로 변환
  colnames(data1) <- names(data)
  data2 <- data.frame(table_name1, data1)
  return(data2)
}))

# mean_rep 테이블 결합 (Wide form으로 열로 저장)
mean_rep <- do.call(rbind, lapply(names(mean_rep_list), function(table_name) {
  data <- mean_rep_list[[table_name]]
  table_name1 <- gsub(pattern = ".RData", "", table_name)
  data1 <- as.data.frame(matrix(data, ncol = 5, byrow = T)) # 값을 행으로 변환
  colnames(data1) <- names(data)
  data2 <- data.frame(table_name1, data1)
  return(data2)
}))

# 테이블 출력
head(mean_nsample)
head(mean_rep)
write.xlsx(mean_nsample, file = "/data4/nekim/matching/Routput_240818/combined_nsample.xlsx", row.names = FALSE)
write.xlsx(mean_rep, file = "/data4/nekim/matching/Routput_240818/combined_rep.xlsx", row.names = FALSE)


