################################
###############################library("devtools")
#setwd("/Volumes/fm6/dNdS/")
#install("dndscv")
library("seqinr")
library("Biostrings")
library("MASS")
library("dndscv")
library(ggplot2)
library(stringr)
library("GenomicRanges")
library(IRanges)
library(mclust)
library(knitr)
library(RColorBrewer)
library(survival)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)

### original mutation file

clin<- read.delim("april_1606_matrix.txt", stringsAsFactors = F)
clin_short<- clin[,c("sampleID","t_CCND1", "t_MMSET","HRD")]

# setwd("~/Desktop/project/SV_project/data/")
# mut<- read.delim("MMRF_CoMMpass_IA11a_IGV_All_Canonical_Variants.mut", stringsAsFactors = F)
# mut2<- mut[,c("sample","chr","start","REF","ALT","TUMOR_ALT_FREQ","GENE")]
# mut2_X<- mut2[mut2$chr=="X",]
# 
# 
# kk<- mut[mut$GENE=="BCL6",]
# kk2<- unique(kk[, c(1:9)])
# 
# length(unique(mut$sample[grep("1_BM",mut$sample)]))
# # [1] 933
# 
# ### upload bad tree strcuture from pyloWGS that we have to remove
# 
# setwd("~/Desktop/project/neutral_evolution//")
# bad_tree<- read.delim("trees_to_exclude_commpass.csv", sep=",", stringsAsFactors = F)
# # dim(bad_tree)
# # [1] 14  2
# 
# ### upload CCF data
# 
# setwd("~/Desktop/project/SV_project/")
# all_tree<- read.delim("CCF_commpass2.txt", sep="\t", stringsAsFactors = F)
# 
# unique(all_tree$chrom) # they eexcluded chromosome X from the analysis
# all_tree<- all_tree[grep("_1_BM",all_tree$sample_ID),]
# # length(unique(all_tree$sample_ID))
# # [1] 802
# 
# all_tree2<- all_tree[!all_tree$sample_ID %in% bad_tree$sample,]
# # length(unique(all_tree2$sample_ID))
# # [1] 788
# all_tree2_select<- all_tree2[,c("sample_ID","chrom","pos","ref","alt")]
# colnames(all_tree2_select) = c("sampleID","chr","pos","ref","mut")
# 
# ### run dnds to correcly annotated each mutation for function, role, effect
# 
# dndsout = dndscv(all_tree2_select)
# mut_ann<- dndsout$annotmuts
# mut_ann$code<- paste(mut_ann$chr, mut_ann$pos, mut_ann$ref, mut_ann$mut, mut_ann$sampleID, sep="_")
# select<- all_tree2[,! colnames(all_tree2) %in% c("chrom","pos","ref","alt")]
# colnames(select)[1]<-"code"
# select$code<- gsub("-","_", select$code)
# 
# # length(unique(select$sample_ID))
# # [1] 788
# 
# ### 76 driver gene on chromosome 1:22
# 
# gne_76<- c("NRAS","FAM46C","LCE1D","ARID1A","CDKN2C","FUBP1","RPL5","NFKB2","ATM", "CCND1","MAML2","PTPN11","DTX1","BCL7A","CDKN1B","KRAS","BHLHE41","ARID2","BTG1",     
#            "RB1", "DIS3","TGDS","TRAF3","NFKBIA","MAX", "ZFP36L1","TCL1A","MAN2C1","IDH2","CREBBP","CYLD","MAF", "NCOR1","TBC1D29","NF1", "TP53","ACTG1","KMT2B" ,   
#            "PRKD2","SF3B1","IDH1","SP140","DNMT3A","DUSP2","SAMHD1","RPRD1B","MAFB","IGLL5","XBP1","EP300","RASA2","RFTN1","PIK3CA","KLHL6","SETD2","TET2","RPS3A" ,   
#            "FGFR3","IRF1","EGR1","PRDM1","HIST1H1E","HIST1H1D","HIST1H2BK","HIST1H1B" , "ABCF1","LTB", "PIM1","IRF4","ZNF292","POT1","BRAF","KMT2C","PABPC1","UBR5","TRAF2" )
# 
# 
# ### upload list of drivers from the merged Maura and Walker list. 
# 
# setwd("~/Desktop/bradley_terry/")
# driver<- read.delim("driver_fra_morgan.txt",sep="\t", stringsAsFactors = F)
# colnames(driver)[1]<-"Gene.Symbol"
# 
# ### find missign drivers on X 
# 
# missing_gene<- driver$Gene.Symbol[! driver$Gene.Symbol %in% gne_76]
# 
# 
# ### merge pylogenetic tree and CCF with driver annotation from dnds
# 
# file1_22<- merge(mut_ann,select, by="code")
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ### 
# ### correct sex/gender 
# ### 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# ## upload gender
# rawdata_all<- read.delim("~/Desktop/bradley_terry/april_1606_matrix.txt", sep="\t", stringsAsFactors = F)
# rawdata_all2<- rawdata_all[,c("sampleID","DEMOG_GENDER")]
# rawdata_all2$sampleID<- paste0(rawdata_all2$sampleID, "_1_BM")
# X_missing<- mut2[mut2$GENE%in% missing_gene,]
# colnames(rawdata_all2)[1]<-"sample"
# X_missing2<- merge(X_missing, rawdata_all2, by="sample")
# 
# ## upload purity
# setwd("~/Desktop/project/NEGATIVE_selection/")
# purity<- read.delim("purity_commpass_a13.txt", sep="\t", stringsAsFactors = F)
# X_missing2_pur<- merge(X_missing2, purity, by="sample") ### CCF is not the mutation CCF but the sample CCF
# 
# 
# ## upload CNV
# setwd("~/Desktop/project/AMoritz/commpass_copy_number_2019/")
# cnv_all<- read.delim("commpass_cnv_new_nov_2019.txt",sep="\t", stringsAsFactors = F)
# length(unique(cnv_all$IDA[grep("1_BM",cnv_all$IDA)]))
# # [1] 812
# 
# cnv_allX<- cnv_all[cnv_all$seqnames=="X",]
# 
# sample_list<- unique(X_missing2_pur$sample)
# all_vaf_X<- list()
# for(i in (1:length(sample_list)))
# {
#   X_missing2_pur_i<- X_missing2_pur[X_missing2_pur$sample == sample_list[i],]
#   cnv_i<- cnv_allX[cnv_allX$IDA == sample_list[i],]
#   gr_driver = with(X_missing2_pur_i, GRanges(chr, IRanges(start=start, end=start)))
#   values(gr_driver) <- DataFrame(REF = X_missing2_pur_i$REF, ntchange = X_missing2_pur_i$ALT, gene= X_missing2_pur_i$GENE, 
#                                  ccf = X_missing2_pur_i$CCF, TUMOR_ALT_FREQ= X_missing2_pur_i$TUMOR_ALT_FREQ, gender= X_missing2_pur_i$DEMOG_GENDER)
#   
#   
#   gr_cnv = with(cnv_i, GRanges(seqnames, IRanges(start=as.numeric(startA), end=as.numeric(endA))))
#   values(gr_cnv) <- DataFrame(IDA=cnv_i$IDA, major=cnv_i$major, minor=cnv_i$minor)
#   
#   range_dri <- merge(as.data.frame(gr_driver),as.data.frame(gr_cnv),by="seqnames",suffixes=c("A","B"))
#   range_dri_2 <- range_dri[with(range_dri, startB <= startA & endB >= endA),]
#   
#   all_vaf_X[[i]]<- range_dri_2
#   
# }
# all_vaf_X2<- do.call("rbind", all_vaf_X)
# all_vaf_X2$mut_ccf<- all_vaf_X2$ccf * all_vaf_X2$TUMOR_ALT_FREQ *all_vaf_X2$major
# all_vaf_X2$mut_ccf[all_vaf_X2$mut_ccf>1]<-1
# plot(density(all_vaf_X2$mut_ccf))
# 
# 
# #### create a loop where CCF of mutations on X chromosome are collapsed to the closest tree cluster 
# 
# sample_list_file<- unique(all_vaf_X2$IDA)
# sample_list_file<- sample_list_file[sample_list_file%in% unique(file$sampleID)]
# x_final<- list()
# for(i in (1:length(sample_list_file)))
# {
#   all_vaf_X2_sam<- all_vaf_X2[all_vaf_X2$IDA== sample_list_file[i],]
#   file_sam<- file[file$sampleID == sample_list_file[i],]
#   CCF_value<- unique(file_sam$CCF)
#   for(j in (1:nrow(all_vaf_X2_sam)))
#   {
#     which(abs(CCF_value - all_vaf_X2_sam$mut_ccf[j]) == min(abs(CCF_value - all_vaf_X2_sam$mut_ccf[j])))
#     all_vaf_X2_sam$mut_ccf[j] <- CCF_value[which(abs(CCF_value - all_vaf_X2_sam$mut_ccf[j]) == min(abs(CCF_value - all_vaf_X2_sam$mut_ccf[j])))]
#   }
#   x_final[[i]]<- all_vaf_X2_sam
# }
# x_final2<- do.call("rbind", x_final)
# x_final2$code<- paste(x_final2$seqnames, x_final2$startA, x_final2$REF, x_final2$ntchange, x_final2$IDA, sep="_")
# x_final2$Clone<- NA
# all_tree2_X<- x_final2[,c("code","Clone","mut_ccf","IDA","seqnames","startA","REF","ntchange")]
# colnames(all_tree2_X)<- colnames(all_tree2)
# 
# x_final2_dnds<- x_final2[,c("IDA","seqnames","startA","REF","ntchange")]
# colnames(x_final2_dnds)<- c("sampleID","chr","pos","ref","mut")
# 
# ### run dnds to get the sme annotation
# 
# dndsout = dndscv(x_final2_dnds)
# mut_ann_X<- dndsout$annotmuts
# mut_ann_X$code<- paste(mut_ann_X$chr, mut_ann_X$pos, mut_ann_X$ref, mut_ann_X$mut, mut_ann_X$sampleID, sep="_")
# select_X<- all_tree2_X[,! colnames(all_tree2_X) %in% c("chrom","pos","ref","alt")]
# colnames(select_X)[1]<-"code"
# select_X$code<- gsub("-","_", select_X$code)
# file_X <- merge(mut_ann_X,select_X, by="code")
# 
# ### create final file
# 
# file_all<- rbind.data.frame(file1_22, file_X) ### this file contains also mutations on X chromosome

# write.table(file_all, "~/Desktop/bradley_terry/rebuttal/file_post_dnds.txt",quote=F, sep="\t")
file_all<- read.delim("~/Desktop/bradley_terry/rebuttal/file_post_dnds.txt",stringsAsFactors = F, sep="\t")

### exclude Synonimuse
file<- file_all[file_all$impact !="Synonymous",]

length(unique(file$sampleID))
# [1] 788
length(unique(file$gene))
# [1] 13351

### upload COSMIC census to define oncosuppressor/oncogene

cosmic<- read.delim("COSMIC_Census.csv",sep=",", stringsAsFactors = F)
cosmic2<- cosmic[,c("Gene.Symbol","Role.in.Cancer")]

### upload List of the 291 AID s discovered in UNG/MSH2 dKO mice

setwd("~/Desktop/bradley_terry/") ### list of 291 AID s in mouse 
aid_jem<- read.delim("aid_targets.txt",sep="\t", stringsAsFactors = F, skip=1)
aid_jem2<- aid_jem[,c("Gene","Ensembl.ID")]

### List of the 18 AID targets mutated in repair-proficient germinal center B cells.

aid_jem_short<- read.delim("aid_targets_short.txt",sep="\t", stringsAsFactors = F)

### List of from Chapuy Nat Med DLBCL Broad paper

setwd("~/Desktop/bradley_terry/nat_med/")
aid_nat_med<- read.delim("41591_2018_16_MOESM6_ESM.txt",sep="\t", stringsAsFactors = F, skip=1)
aid_nat_med_sig<- aid_nat_med[aid_nat_med$qval<0.05,]
aid_nat_med2<- aid_nat_med_sig[,c("gene","qval")]
colnames(aid_nat_med2)<- colnames(aid_jem2)
nrow(aid_nat_med2)
# [1] 33

### we selected the more comprehensive list

aid<- unique(rbind.data.frame(aid_jem2, aid_nat_med2))

aid_short<- unique(c(aid_jem_short$Gene, aid_nat_med2$Gene))

### function to convert all gene in names in capital letter

target<- data.frame(lapply(aid, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))
colnames(target)[1]<-"Gene.Symbol"
target$aid<- "yes"


### create final file with COSMIC TSG/Oncogene, AID annotation
# driver is the list of 81 drivers 
require(plyr)
driver$Gene.Symbol<- gsub("[*]","", driver$Gene.Symbol)
driver_fun<- join( driver, cosmic2,by="Gene.Symbol") ### join list of 81 driver mutations with COSMIC annotation 
driver_fun2<- driver_fun[,c("Gene.Symbol","Role")]
functio<- unique(join(driver_fun2, target[,-2], by="Gene.Symbol")) ### join with AID taregt information

### select nonsyno. mut in driver genes. 

final_mut<- file[file$gene %in% driver[,1],]


### check numbers

length(unique(final_mut$sampleID))
# [1] 700
length(unique(final_mut$gene))
# [1] 81

#### select nonsyno. mut not involving driver genes. 

final_no_driver<- file[!file$gene %in% driver[,1],]

###################################################
###
### Bradley Terry Model analysis
###
#####################################################

set1 <- brewer.pal(8, "Set1")
dark2 <- brewer.pal(8, "Dark2")

### cerate data.frame for bradley terry model following original Elli NEJM apporach

elli3<- final_mut[,c("sampleID","CCF","Clone","gene","CCF")] ### keep the columns as the original data frame to avoid issues
colnames(elli3)[1:5]<-c("Sample","PM.Tum","DP.Tum","Gene","BT_VAF")
colnames(elli3)[5]<-"VAF_adj"
elli3$Sample<-as.character(as.character(elli3$Sample))
elli3$Gene<-as.character(as.character(elli3$Gene))

### some passages were included by elli to correct targeted sequencing possibile bias.

sam_cyto_group<- list()
for(i in (2:4))
{
  sam_cyto_group[[colnames(clin_short[i])]]<- unique(clin_short[clin_short[,i] == 1,]$sampleID)
}

for(ww in (1:3))
{
  mutationData<- elli3[elli3$Sample %in% paste0(sam_cyto_group[[ww]],"_1_BM"),]
print(length(unique(mutationData$Sample)))
  }
  
for(ww in (1:3))
{
mutationData<- elli3[elli3$Sample %in% paste0(sam_cyto_group[[ww]],"_1_BM"),]

mutationData<- elli3[elli3$Sample %in% c(paste0(sam_cyto_group[[2]],"_1_BM"),paste0(sam_cyto_group[[1]],"_1_BM")),]

mutationData<- mutationData[order(mutationData$Sample),]
mutationTable <- (table(mutationData[c("Sample","Gene")]) > 0)+0

mutationData$DP.Tum <- as.numeric(as.character(mutationData$DP.Tum))
vafCn <- as.numeric(as.character(mutationData$VAF_adj))
vaf <- mutationData$PM.Tum 
mutationData$Sample <- factor(as.character(mutationData$Sample))
vafCn <- as.numeric(as.character(mutationData$VAF_adj))

list_seg<- unique(mutationData$Gene)
clinicalData<-as.data.frame.matrix(table(elli3$Sample,elli3$Gene))

### pairwise precedence in patients with more than 1 driver mutations

precedence <- matrix(0, ncol=ncol(mutationTable), nrow=ncol(mutationTable), dimnames=list(colnames(mutationTable), colnames(mutationTable)))
plist <- list()

ix=mutationData$Gene %in% colnames(precedence)
list_name<- as.character(unique(mutationData$Sample))
for(s in list_name){
  l <- list()
  for(i in which(mutationData$Sample==s & ix))
    for(j in which(mutationData$Sample==s & ix)){
      if(i!=j){
        if(vafCn[i] > vafCn[j]){ 
          #Pidgeonhole
          precedence[as.character(mutationData$Gene[i]),as.character(mutationData$Gene[j])] <- precedence[as.character(mutationData$Gene[i]),as.character(mutationData$Gene[j])] + 1
          l <- c(l, list(c(as.character(mutationData$Gene[i]),as.character(mutationData$Gene[j]))))       
        }
      }
    }
  plist[[s]] <- l
}


sample_list<- unique(elli3$sample)

### 2 functions for bradley terry model, to generate IC

makeDesign <- function(I) {
  w <- which(lower.tri(I), arr.ind=TRUE)
  x <- matrix(0, nrow(w), nrow(I))
  for(i in 1:nrow(w)){
    x[i,w[i,1]] <- 1
    x[i,w[i,2]] <- -1
  }
  return(x)
}

btModel <- function(I){
  y <- cbind(I[lower.tri(I)], t(I)[lower.tri(I)])
  x <- makeDesign(I = I)
  glm.fit(x=x[,-1],y=y, family=binomial())
}

precedence2<- precedence

setwd("~/Desktop/bradley_terry/")
# pdf("bradley_terry_integrated_all_mol_time_order_new_copy-number.pdf", width=8, heigh=8)
par(mfrow=c(1,1))
par(mar=c(3,3,3,0),xpd=T )
nCasesGene <- table(factor(unlist(sapply(plist, function(x) unique(unlist(x)))), levels=colnames(precedence)))
# w <- which(nCasesGene > 7)
w<- which(nCasesGene > 1)
fit <- btModel(precedence[w,w]+.01) ## Warning: non-integer counts in a binomial glm!
c <- c(0,coef(fit))
names(c) <- colnames(precedence)[w]
o <- rank(c)
v <- pmin(2,sqrt(c(0,diag(chol2inv(fit$qr$qr)))))
l <- names(c)
m <- paste("n=",nCasesGene[w], sep="")

functio<- join(driver_fun2, target, by="Gene.Symbol") ### import AID and COSMIC info
functio<- unique(functio[,c(1,2,4)])
functio$aid[is.na(functio$aid)]<-"brown2"
functio$aid[functio$aid=="yes"]<-"dodgerblue"
functio$Role[functio$Role=="TSG"]<-"forestgreen"
functio$Role[functio$Role==""]<-"forestgreen"
functio$Role[functio$Role=="Oncogene"]<-"darkorchid4"
functio$Role[functio$Role=="Unknown"]<-"grey20"
functio<- functio[order(functio$Gene.Symbol),]
functio2<- functio[functio$Gene.Symbol%in%names(c),]

setwd("~/Desktop/bradley_terry/rebuttal/figure_rebuttal/")
pdf("IGH_bradley_terry_aid_annotation.pdf", width=8, heigh=10)
# pdf(paste(names(sam_cyto_group[ww]),"bradley_terry_aid_annotation.pdf"), width=8, heigh=10)
par(mfrow=c(1,1))
par(mar=c(3,3,3,0),xpd=T )
plot(-c, o, xlab="Relative time", yaxt="n", pch=19, col=functio$aid[match(names(c),functio$Gene.Symbol )], ylab="", xlim=c(-8,10), bty="n")
segments(-c-v, o,-c+v,o, col="grey")
par(new=TRUE)
plot(-c, o, xlab="", yaxt="n",xaxt="n", pch=19, col=functio$aid[match(names(c),functio$Gene.Symbol )],
     ylab="", xlim=c(-8,10), bty="n")
text(-c-v ,o,gsub("_"," ",l), font=3, pos=2, cex=0.8,  col=functio$Role[match(names(c),functio$Gene.Symbol )])
text(-c+v ,o,m, font=1, pos=4, cex=1)
dev.off()
}
