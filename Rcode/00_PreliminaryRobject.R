rm(list=ls())
.libPaths("/data/pkg36")
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/DealFunc.R")
source("/home2/nekim/scratch/matching/Rcode/Rcode_DegreeTotal/SourceCode_Binary_Scene1_2.R") # 
load(file="/home2/nekim/scratch/matching/Robject_old/seed.RData")

# metagenome 코드 불러오기
physeq<-as.data.frame(fread("/home2/nekim/scratch/NOTOfull/output/otu_table_mc2_w.txt"))
tree <- read.tree(file="/home2/nekim/scratch/NOTOfull/output/rep_set_mc2.tre")
colnames(physeq)
otu <- otu_table(physeq[,2:(ncol(physeq)-1)], taxa_are_rows = T)
rownames(otu) <- physeq[,1]
tax <- tax_table(as.matrix(data.frame(tx=physeq[,ncol(physeq)]) %>% separate(tx, c("Domain","Phylum","Class","Order","Family","Genus","Species"), ";")))
rownames(tax) <- physeq[,1]
phy <- phyloseq(otu_table(otu), tax_table(tax))
any(taxa_sums(phy)==0)
ntaxa(phy)

nsamp=1000
set.seed(1)
toysam <- c(rep(TRUE, 1000), rep(FALSE, nsamples(phy)-1000))
phyy <- prune_samples(toysam, phy)
otutable <- filtOTU(otu_table(phyy))
indicator=indc(otutable, omics="Metagenomics")
lbsize <- apply(phyy@otu_table,2,sum)

toysam200 <- c(rep(TRUE, 200), rep(FALSE, nsamples(phy)-200)) 
phyy200 <- prune_samples(toysam200, phy)
otutable200 <- filtOTU(otu_table(phyy200))
indicator200=indc(otutable200, omics="Metagenomics")

toysam500 <- c(rep(TRUE, 500), rep(FALSE, nsamples(phy)-500)) 
phyy500 <- prune_samples(toysam500, phy)
otutable500 <- filtOTU(otu_table(phyy500))
indicator500=indc(otutable500, omics="Metagenomics")

otulist <- list(otutable, otutable500, otutable200, lbsize)
indiclist <- list(indicator, indicator500, indicator200)

save(otulist, file="/home2/nekim/scratch/matching/Robject/otulist.RData")
save(indiclist, file="/home2/nekim/scratch/matching/Robject/indiclist.RData")


psdiff <- unlist(lapply(data.frame(t(otutable)), psdiffun, k=indicator))
psdiff500 <- unlist(lapply(data.frame(t(otutable500)), psdiffun, k=indicator500))
psdiff200 <- unlist(lapply(data.frame(t(otutable200)), psdiffun, k=indicator200))

psdifflist <- list(psdiff, psdiff500, psdiff200)
save(psdifflist, file="/home2/nekim/scratch/matching/Robject/psdifflist.RData")
