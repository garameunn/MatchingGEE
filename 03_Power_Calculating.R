
rm(list=ls())
options(digits = 3)
.libPaths("/data/pkg36")
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

## 사실상 scenario 1,2,3 만 다르게 입력하면 똑같이 돌릴 수 있음.

#####################################
## 챕터 : 1 (Continuous)
## 시나리오: 1,2,3
## 통계량 : Power
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : SourceCode_Continuous_Scene123.R
#####################################
#### nsamp, scenario 변경
nsamp = 200 # 200, 500, 1000

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]] ## 이건 아니지!!
lbs = otulist[[4]][1:nsamp]

b = c(0.001,0.002,0.005, 0.01)
k <- c(1,8,61,77)

scenario="S2" # S1, S2, S3

#### rate, rep, nonlin 변경
set.seed(1)
rep=500
dimm=dim(otutable)[1]; omics="Metagenomics"
datatable=otutable;rate=0.2;nonlin="inv" # comp
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Continuous_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 
limtable <- removeBatchEffect(trsftable, batch=bat.dat1$Batch)

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    if (nonlin == "inv"){
      tablefun <- get("tablemaking_powerC")
      tablefun1 <- get("tablemaking_power_woC")
    } else if (nonlin == "comp") {
      tablefun <- get("tablemaking_power.comp")
      tablefun1 <- get("tablemaking_power.comp_wo")
    }
    
    
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
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- limm.sig <- limm.coef <- c()
    
    
    for (i in 1:length(b)){
      
      tble <- data.frame(tble0[[1]][[i]])
      tble0.1 <- tble0_.1[[1]][[i]]
      tble0.2 <- tble0_.2[[1]][[i]]
      tble0.3 <- tble0_.3[[1]][[i]]
      tble0.05 <- tble0_.05[[1]][[i]]
      tble_wo <- tble0_wo[[1]][[i]]
      tble0.1_wo <- tble0_.1_wo[[1]][[i]]
      tble0.2_wo <- tble0_.2_wo[[1]][[i]]
      tble0.3_wo <- tble0_.3_wo[[1]][[i]]
      tble0.05_wo <- tble0_.05_wo[[1]][[i]]
      
      tble2 <- data.frame(tble0[[2]][[i]]) # entire data
      tble2$limma <- unlist(limtable[x,])
      tble2$cmbt <- unlist(cmbatable[x,])
      
      
      #######################
      ###### with rep #######
      #######################
      # no caliper
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig <- c(r.sig,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff <- c(coeff,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig2 <- c(r.sig2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff2 <- c(coeff2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      rm(geemfit);rm(geemfit2); 
      
      # 0.2 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.2 <- c(r.sig0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.2 <- c(coeff0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.2 <- c(r.sig20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.2 <- c(coeff20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      rm(geemfit);rm(geemfit2); 
      
      # 0.1 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.1 <- c(r.sig0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.1 <- c(coeff0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.1 <- c(r.sig20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.1 <- c(coeff20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      rm(geemfit);rm(geemfit2);
      
      # 0.3 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.3 <- c(r.sig0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.3 <- c(coeff0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.3 <- c(r.sig20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.3 <- c(coeff20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.05 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.05 <- c(r.sig0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.05 <- c(coeff0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.05 <- c(r.sig20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.05 <- c(coeff20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      
      #########################
      ###### without rep ######
      #########################
      rm(geemfit);rm(geemfit2); 
      # no caliper
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig_w <- c(r.sig_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff_w <- c(coeff_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig2_w <- c(r.sig2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff2_w <- c(coeff2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.2 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.2_w <- c(r.sig0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.2_w <- c(coeff0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.2_w <- c(r.sig20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.2_w <- c(coeff20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.1 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.1_w <- c(r.sig0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.1_w <- c(coeff0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.1_w <- c(r.sig20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.1_w <- c(coeff20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      # 0.3 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.3_w <- c(r.sig0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.3_w <- c(coeff0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.3_w <- c(r.sig20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.3_w <- c(coeff20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.05 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="gaussian", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.05_w <- c(r.sig0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.05_w <- c(coeff0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="gaussian", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.05_w <- c(r.sig20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.05_w <- c(coeff20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      rm(geemfit);rm(geemfit2);
      
      
      logit <- lm(paste("case~x1+PS","marker",sep="+"), data=tble2)
      logit.sig <- ifelse(abs(summary(logit)$coefficients[4,1])<100,summary(logit)$coefficients[4,4],1)
      logit.coef <- ifelse(abs(summary(logit)$coefficients[4,1])<100,coef(summary(logit))[4,1],NA)
      rm(logit)
      
      logit1 <- lm(paste("case~x1","marker",sep="+"), data=tble2)
      logit.sig1 <- ifelse(abs(summary(logit1)$coefficients[3,1])<100,summary(logit1)$coefficients[3,4],1)
      logit.coef1 <- ifelse(abs(summary(logit1)$coefficients[3,1])<100,coef(summary(logit1))[3,1],NA)
      rm(logit1)
      
      
      combat <- lm(paste("case~x1+hetero1+hetero2","cmbt",sep="+"), data=tble2) 
      combat.sig <- ifelse(abs(summary(combat)$coefficients[5,1])<100,summary(combat)$coefficients[5,4],1)
      combat.coef <- ifelse(abs(summary(combat)$coefficients[5,1])<100,coef(summary(combat))[5,1],NA)
      
      limm <- lm(paste("case~x1+hetero1+hetero2","limma",sep="+"), data=tble2) 
      limm.sig <- ifelse(abs(summary(limm)$coefficients[5,1])<100,summary(limm)$coefficients[5,4],1)
      limm.coef <- ifelse(abs(summary(limm)$coefficients[5,1])<100,coef(summary(limm))[5,1],NA)
      
    }
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, limm.sig,
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #24
                      logit.coef1, logit.coef, combat.coef, limm.coef,
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w)) #24   
    
    
  },mc.cores=8))
},mc.cores = 6)
head(pvalues[[1]])

save(pvalues, file=paste0("/data4/nekim/matching/Robject/Power_conti/NSAMP",nsamp,"_REP",rep,"_",nonlin,"_",scenario,".RData")) 












#####################################
## 챕터 : 2 (Binary)
## 시나리오: 1
## 통계량 : Power
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : tmp_power.R
#####################################

#### nsamp 변경
nsamp = 200

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = c(0.001,0.002,0.005, 0.01)

scenario="S1"

#### rate, rep, nonlin 변경
set.seed(1)
rep=500
dimm=dim(otutable)[1]; omics="Metagenomics"
datatable=otutable;rate=0.5;nonlin="comp" # comp
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 
limtable <- removeBatchEffect(trsftable, batch=bat.dat1$Batch)

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    
    tablefun <- get("tablemaking_power")
    tablefun1 <- get("tablemaking_power_wo")
    
    
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
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- limm.sig <- limm.coef <- c()
    
    
    for (i in 1:length(b)){
      
      tble <- data.frame(tble0[[1]][[i]])
      tble0.1 <- tble0_.1[[1]][[i]]
      tble0.2 <- tble0_.2[[1]][[i]]
      tble0.3 <- tble0_.3[[1]][[i]]
      tble0.05 <- tble0_.05[[1]][[i]]
      tble_wo <- tble0_wo[[1]][[i]]
      tble0.1_wo <- tble0_.1_wo[[1]][[i]]
      tble0.2_wo <- tble0_.2_wo[[1]][[i]]
      tble0.3_wo <- tble0_.3_wo[[1]][[i]]
      tble0.05_wo <- tble0_.05_wo[[1]][[i]]
      
      tble2 <- data.frame(tble0[[2]][[i]]) # entire data
      tble2$limma <- unlist(limtable[x,])
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
      
      limm <- glm(paste("case~x1+hetero1+hetero2","limma",sep="+"), data=tble2, family = "binomial") 
      limm.sig <- c(limm.sig,ifelse(limm$converged,coef(summary(limm))[5,4],1))
      limm.coef <- c(limm.coef,ifelse(limm$converged,coef(summary(limm))[5,1],NA))

    }
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, limm.sig,
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #24
                      logit.coef1, logit.coef, combat.coef, limm.coef,
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w)) #24   
    
    
  },mc.cores=8))
},mc.cores = 6)
head(pvalues[[1]])

save(pvalues, file=paste0("/data4/nekim/matching/Robject/Power/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 
load("/data4/nekim/matching/Robject/Power/NSAMP1000_REP500_RATE0.2_inv.RData")




#####################################
## 챕터 : 2 (Binary)
## 시나리오: 2
## 통계량 : Power
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : OmicGee_ver14.R
#####################################

#### nsamp 변경
nsamp = 500

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = c(0.001,0.002,0.005, 0.01)

scenario="S2"


#### rate, rep, nonlin 변경
set.seed(1)
rep=500
dimm=dim(otutable)[1]
datatable=otutable;rate=0.5;nonlin=NULL
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 
limtable <- removeBatchEffect(trsftable, batch=bat.dat1$Batch)

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    if (nonlin == "inv"){
      tablefun <- get("tablemaking_power")
      tablefun1 <- get("tablemaking_power_wo")
    } else if (nonlin == "comp") {
      tablefun <- get("tablemaking_power.comp")
      tablefun1 <- get("tablemaking_power.comp_wo")
    }
    
    
    tble0 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, nonlin=nonlin)  
    tble0_.05 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, nonlin=nonlin)
    tble0_.1 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, nonlin=nonlin)
    tble0_.2 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, nonlin=nonlin)
    tble0_.3 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, nonlin=nonlin)
    
    tble0_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, nonlin=nonlin) 
    tble0_.05_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, nonlin=nonlin)  
    tble0_.1_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, nonlin=nonlin)
    tble0_.2_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, nonlin=nonlin)  
    tble0_.3_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, nonlin=nonlin)
    
    
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
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- limm.sig <- limm.coef <- c()
    
    
    for (i in 1:length(b)){
      
      tble <- data.frame(tble0[[1]][[i]])
      tble0.1 <- tble0_.1[[1]][[i]]
      tble0.2 <- tble0_.2[[1]][[i]]
      tble0.3 <- tble0_.3[[1]][[i]]
      tble0.05 <- tble0_.05[[1]][[i]]
      tble_wo <- tble0_wo[[1]][[i]]
      tble0.1_wo <- tble0_.1_wo[[1]][[i]]
      tble0.2_wo <- tble0_.2_wo[[1]][[i]]
      tble0.3_wo <- tble0_.3_wo[[1]][[i]]
      tble0.05_wo <- tble0_.05_wo[[1]][[i]]
      
      tble2 <- data.frame(tble0[[2]][[i]]) # entire data
      tble2$limma <- unlist(limtable[x,])
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
      
      limm <- glm(paste("case~x1+hetero1+hetero2","limma",sep="+"), data=tble2, family = "binomial") 
      limm.sig <- c(limm.sig,ifelse(limm$converged,coef(summary(limm))[5,4],1))
      limm.coef <- c(limm.coef,ifelse(limm$converged,coef(summary(limm))[5,1],NA))
      
    }
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, limm.sig,
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #24
                      logit.coef1, logit.coef, combat.coef, limm.coef,
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w)) #24   
    
    
  },mc.cores=8))
},mc.cores = 6)
head(pvalues[[1]])

save(pvalues, file=paste0("/data4/nekim/matching/Robject/Power/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData"))  # 원래 scenario 2,3 에서도 nonlin 이 적용 안되더라도 inv 로 넣고 했었는데, NULL 로 수정한 상태임




#####################################
## 챕터 : 2 (Binary)
## 시나리오: 3
## 통계량 : Power
## 주제 : 계산 (계산/정리 중 하나)
## 출처 : OmicGee_ver14.R
#####################################

#### nsamp 변경
nsamp = 200

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = c(0.001,0.002,0.005, 0.01)

scenario="S3"


# k 설정
k <- c(1,8,61,77)



#### rate, rep, nonlin 변경
set.seed(1)
rep=500
dimm=dim(otutable)[1]
datatable=otutable;rate=0.5;nonlin=NULL
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene123.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 
limtable <- removeBatchEffect(trsftable, batch=bat.dat1$Batch)

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    if (nonlin == "inv"){
      tablefun <- get("tablemaking_power")
      tablefun1 <- get("tablemaking_power_wo")
    } else if (nonlin == "comp") {
      tablefun <- get("tablemaking_power.comp")
      tablefun1 <- get("tablemaking_power.comp_wo")
    }
    
    
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
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- limm.sig <- limm.coef <- c()
    
    
    for (i in 1:length(b)){
      
      tble <- data.frame(tble0[[1]][[i]])
      tble0.1 <- tble0_.1[[1]][[i]]
      tble0.2 <- tble0_.2[[1]][[i]]
      tble0.3 <- tble0_.3[[1]][[i]]
      tble0.05 <- tble0_.05[[1]][[i]]
      tble_wo <- tble0_wo[[1]][[i]]
      tble0.1_wo <- tble0_.1_wo[[1]][[i]]
      tble0.2_wo <- tble0_.2_wo[[1]][[i]]
      tble0.3_wo <- tble0_.3_wo[[1]][[i]]
      tble0.05_wo <- tble0_.05_wo[[1]][[i]]
      
      tble2 <- data.frame(tble0[[2]][[i]]) # entire data
      tble2$limma <- unlist(limtable[x,])
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
      
      limm <- glm(paste("case~x1+hetero1+hetero2","limma",sep="+"), data=tble2, family = "binomial") 
      limm.sig <- c(limm.sig,ifelse(limm$converged,coef(summary(limm))[5,4],1))
      limm.coef <- c(limm.coef,ifelse(limm$converged,coef(summary(limm))[5,1],NA))
      
    }
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, limm.sig,
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #24
                      logit.coef1, logit.coef, combat.coef, limm.coef,
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w)) #24   
    
    
  },mc.cores=8))
},mc.cores = 6)
head(pvalues[[1]])

save(pvalues, file=paste0("/data4/nekim/matching/Robject/Power/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 


















#####################################
#####################################
## 챕터 : 1,2 
## 시나리오: 1,2,3
## 통계량 : Power
## 주제 : 정리 (계산/정리 중 하나)
## 출처 : tmp_power.R
#####################################
source("/data4/nekim/matching/Rcode/Rcode_DegreeTotal/DealFunc.R")
load(file="/data4/nekim/matching/Robject/EffectSizeList.RData") # b2list

Binary=TRUE
for (i in c(1,2,3)){
  for (j in c("S1","S2","S3")){
    if (j=="S1"){
      for (nonlin in c("comp","inv")){
        for (Rate in c(0.2,0.3,0.5)){
          tablesfun(norder=i, scene=j, bin=Binary, prevalence=Rate, nonlin=nonlin)
        }
      }
    } else {
      for (Rate in c(0.2,0.3,0.5)){
        tablesfun(norder=i, scene=j, bin=Binary, prevalence=Rate, nonlin="inv")
      }
    }

  }
}

# 
# 
# 
# nsamp<-500;scenario<-"S1";Rate="0.3";rep=500;nonlin="inv"
# 
# load(file=paste0("/data4/nekim/matching/Robject/Power/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 
# 
# 
# 
# 
# 
# 
# 
# ###############
# ## 정리 코드 ##
# ###############
# # load(file="/data4/nekim/matching/Robject/power_5000_nominalB_summary.RData") # nblist=list(nb2_c1s1, nb2_top20)
# load(file="/data4/nekim/matching/Robject/power_5000_nominalBNEW_summary.RData") # nblist=list(nb2_c1s1, nb2_top20)
# 
# rate=0.3
# 
# # robust
# load(file=paste0("/data4/nekim/matching/Robject/Power/rawp_500_Y",rate,"_robust.RData")) # rep=500, dimm=80
# pw <- pvalues
# 
# # naive
# load(file=paste0("/data4/nekim/matching/Robject/Power/rawp_500_Y",rate,".RData")) # rep=500, dimm=80
# load(file=paste0("/data4/nekim/matching/Robject/Power/rawp_500_Y",rate,"NEW.RData")) # rep=500, dimm=80
# 
# pw <- pvalues
# 
# 
# 
# powerfold <- function(p){
#   # browser()
#   p0001 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[1,]))), byrow=T, ncol=13)
#   p0005 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[2,]))), byrow=T, ncol=13)
#   p001 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[3,]))), byrow=T, ncol=13)
#   p002 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[4,]))), byrow=T, ncol=13)
#   return(list(p0001,p0005,p001,p002))
# }
# powerfold1 <- function(x){
#   t005 <- sum(ifelse(is.na(x),1,x)<0.05)/length(x)
# }
# 
# pvalues_p <- lapply(pw, function(x) lapply(x, function(xx) xx[,1:13])) # 1:13=sig, 14:26=beta
# pvalues1 <- powerfold(pvalues_p)
# powers <- do.call(rbind,lapply(pvalues1, function(xx) apply(xx, 2, powerfold1)))
# 
# rownames(powers) <-c("var.0001","var.0005","var.001","var.002") # row=marker variance % 가 0.0001~0.002 일 때 power
# colnames(powers) <- c("covariate_adjusted","wrep_nocal_logit","wrep_nocal_ex","wrep_cal0.2_logit","wrep_cal0.2_ex","wrep_cal0.1_logit","wrep_cal0.1_ex",
#                       "worep_nocal_logit","worep_nocal_ex","worep_cal0.2_logit","worep_cal0.2_ex","worep_cal0.1_logit","worep_cal0.1_ex")
# 
# conf.int <- do.call(rbind,lapply(pvalues1, function(y) apply(y, 2, FUN=CIfun, level=0.05)))
# rownames(conf.int) <- paste0("CI_",c(rep(c("v.0001","v.0005","v.001","v.002"),each=2)))
# powers1 <- rbind(powers, conf.int)
# t(powers1)
# powers2 <-powers1
# 
# 
# 
# powerBfold <- function(p, diffbig){
#   # browser()
#   p005 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[1,]))), byrow=T, ncol=13)
#   p01 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[2,]))), byrow=T, ncol=13)
#   p02 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[3,]))), byrow=T, ncol=13)
#   p04 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[4,]))), byrow=T, ncol=13)
#   return(list(p005,p01,p02,p04))
# }
# betas_p <- lapply(pw, function(x) lapply(x, function(xx) xx[,14:26])) # 1:13=sig, 14:26=beta
# 
# betas1 <- powerBfold(betas_p)
# betamean <- matrix(unlist(lapply(betas1, function(x) { # 그냥 평균
#   d <- apply(x,2,mean, na.rm=T)
#   return(d)
# })),ncol=13, byrow=T)
# colnames(betamean) <- c("covariate_adjusted","wrep_nocal_logit","wrep_nocal_ex","wrep_cal0.2_logit","wrep_cal0.2_ex","wrep_cal0.1_logit","wrep_cal0.1_ex",
#                         "worep_nocal_logit","worep_nocal_ex","worep_cal0.2_logit","worep_cal0.2_ex","worep_cal0.1_logit","worep_cal0.1_ex")
# betamean # 모든 방법에서 OTU 평균적으로 추정된 beta 추정치
# betamean1 <- betamean
# 
# 
# # Top 20
# 
# pickp <- function(x,otus) {
#   xx <- x[otus]
#   return(xx)
# }
# 
# pvalues_pp<-lapply(pvalues_p, pickp, otus=order(psdiff, decreasing = T)[1:20])
# 
# pvalues1 <- powerfold(pvalues_pp)
# powers <- do.call(rbind,lapply(pvalues1, function(xx) apply(xx, 2, powerfold1)))
# 
# rownames(powers) <-c("var.0001","var.0005","var.001","var.002") # row=marker variance % 가 0.0001~0.002 일 때 power
# colnames(powers) <- c("covariate_adjusted","wrep_nocal_logit","wrep_nocal_ex","wrep_cal0.2_logit","wrep_cal0.2_ex","wrep_cal0.1_logit","wrep_cal0.1_ex",
#                       "worep_nocal_logit","worep_nocal_ex","worep_cal0.2_logit","worep_cal0.2_ex","worep_cal0.1_logit","worep_cal0.1_ex")
# 
# conf.int <- do.call(rbind,lapply(pvalues1, function(y) apply(y, 2, FUN=CIfun, level=0.05)))
# rownames(conf.int) <- paste0("CI_",c(rep(c("v.0001","v.0005","v.001","v.002"),each=2)))
# powers1 <- rbind(powers, conf.int)
# t(powers1)
# 
# 
# 
# betas_pp<-lapply(betas_p, pickp, otus=order(psdiff, decreasing = T)[1:20])
# 
# betas1 <- powerBfold(betas_pp)
# betamean <- matrix(unlist(lapply(betas1, function(x) { # 그냥 평균
#   d <- apply(x,2,mean, na.rm=T)
#   return(d)
# })),ncol=13, byrow=T)
# colnames(betamean) <- c("covariate_adjusted","wrep_nocal_logit","wrep_nocal_ex","wrep_cal0.2_logit","wrep_cal0.2_ex","wrep_cal0.1_logit","wrep_cal0.1_ex",
#                         "worep_nocal_logit","worep_nocal_ex","worep_cal0.2_logit","worep_cal0.2_ex","worep_cal0.1_logit","worep_cal0.1_ex")
# betamean # 모든 방법에서 OTU 평균적으로 추정된 beta 추정치
# 
# 
# 
# 
# t(betamean1) #  total
# t(betamean) # top20
# 
# betamean1
# nblist[[1]]
# t(nblist[[1]]) # total
# t(nblist[[2]]) # top 20
# 
# 
# betamean1[,1]
# 
# dtfun <- function(dt, ref){
#   col1 <- (dt[,1]-ref[,1])/ref[,1]; col2 <- (dt[,2]-ref[,2])/ref[,2]; col3 <- (dt[,3]-ref[,2])/ref[,2]; 
#   col4 <- (dt[,4]-ref[,3])/ref[,3]; col5 <- (dt[,5]-ref[,3])/ref[,3]; 
#   col6 <- (dt[,6]-ref[,4])/ref[,4]; col7 <- (dt[,7]-ref[,4])/ref[,4]; 
#   col8 <- (dt[,8]-ref[,5])/ref[,5]; col9 <- (dt[,9]-ref[,5])/ref[,5]; 
#   col10 <- (dt[,10]-ref[,6])/ref[,6]; col11 <- (dt[,11]-ref[,6])/ref[,6]; 
#   col12 <- (dt[,11]-ref[,7])/ref[,7]; col13 <- (dt[,13]-ref[,7])/ref[,7]; 
#   dtt <- data.frame(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13)
#   return(dtt)
# }
# relbias1 <- dtfun(betamean1, nblist[[1]]) # total
# relbias <- dtfun(betamean, nblist[[2]]) # top20
# colnames(relbias)<-colnames(relbias1) <- c("covariate_adjusted","wrep_nocal_logit","wrep_nocal_ex","wrep_cal0.2_logit","wrep_cal0.2_ex","wrep_cal0.1_logit","wrep_cal0.1_ex",
#                                            "worep_nocal_logit","worep_nocal_ex","worep_cal0.2_logit","worep_cal0.2_ex","worep_cal0.1_logit","worep_cal0.1_ex")
# t(relbias1) # total
# t(relbias) # top20
# 
# 
# # relbias1 <- dtfun(exp(betamean1), exp(nblist[[1]])) # total
# # relbias <- dtfun(exp(betamean), exp(nblist[[2]])) # top20
# # colnames(relbias)<-colnames(relbias1) <- c("covariate_adjusted","wrep_nocal_logit","wrep_nocal_ex","wrep_cal0.2_logit","wrep_cal0.2_ex","wrep_cal0.1_logit","wrep_cal0.1_ex",
# #                                            "worep_nocal_logit","worep_nocal_ex","worep_cal0.2_logit","worep_cal0.2_ex","worep_cal0.1_logit","worep_cal0.1_ex")
# # t(relbias1)*100 # total
# # t(relbias)*100 # top20
# 
# 
# 
# # plot
# powers1 # top20
# powers2 # total
# 
# 
# pd <- powers1[1:4,c(1,3,5,9,11)]
# rownames(pd) <- paste0("v",c("01","02","05", "1"))
# colnames(pd) <- c("M1","M2","M3","M4","M5") # "covariate_adjusted" "wrep_nocal_ex"      "wrep_cal0.2_ex"     "worep_nocal_ex"     "worep_cal0.2_ex"  
# pdm <- melt(pd)
# png(file=paste0("/data4/nekim/matching/Robject/Power/T1_scene1_Yper",rate,"_top20.png"), width=300,height=450)
# p<-ggplot(pdm, aes(x=X1, y=value, group=X2)) +
#   geom_line(aes(color=X2, linetype=X2))+
#   geom_point(aes(color=X2))+ theme_bw() + ylim(c(0.05,0.6)) +
#   xlab("Effect size")+ylab("Power")+
#   labs(color = "Model", linetype="Model", title = paste0("[Y=1] proportion: ",rate*100,"%"))
# print(p)
# dev.off()
# 
# 
# 
# 
# relbias1 # total
# relbias # top20
# rownames(relbias) <- rownames(relbias1) <- paste0("v",c("01","02","05", "1"))
# 
# 
# # rd <- relbias1[1:4,c(1,3,5,9,11)]
# rd <- relbias1[1:4,c(1,3,2,5,4,9,8,11,10)]
# rownames(rd) <- paste0("v",c("01","02","05", "1"))
# # colnames(rd) <- c("M1","M2","M3","M4","M5") # "covariate_adjusted" "wrep_nocal_ex"      "wrep_cal0.2_ex"     "worep_nocal_ex"     "worep_cal0.2_ex"  
# 
# colnames(rd) <- c("M1","M2","M2-1","M3","M3-1","M4","M4-1","M5","M5-1") # covariate_adjusted" "wrep_nocal_ex"      "wrep_nocal_logit"   "wrep_cal0.2_ex"     "wrep_cal0.2_logit"  "worep_nocal_ex"     "worep_nocal_logit"  "worep_cal0.2_ex"    "worep_cal0.2_logit"
# rd<-as.matrix(rd)
# rdm <- melt(rd)
# png(file=paste0("/data4/nekim/matching/Robject/Power/T1_scene1_Yper",rate,"_total_RB.png"), width=300,height=450)
# pp<-ggplot(rdm, aes(x=X1, y=value, group=X2)) +
#   geom_line(aes(color=X2, linetype=X2))+
#   geom_point(aes(color=X2))+ theme_bw() +ylim(c(-1.5,1.5)) +
#   xlab("Nominal Effect size")+ylab("Relative Bias")+
#   labs(color = "Model", linetype="Model", title = paste0("[Y=1] proportion: ",rate*100,"%"))
# print(pp)
# dev.off()
# 

