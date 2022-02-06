#### plot Mol_time comparing 

### CNV overlaps with GISTIC peaks
library(stringr)




icgc_grch37<-read.table(file = "PATH/TO/sample_normal_contamination.txt", header = TRUE)
gambl_grch37_TN<-as.data.frame(str_split_fixed(gambl_grch37$SAMPLE, "--", 2))
colnames(gambl_grch37_TN)<-c("Tumour","Normal")
gambl_grch37<-cbind(gambl_grch37,gambl_grch37_TN)


all<-rbind(icgc_grch37,gambl_grch37)

## metadata file for the study
#fl_meta_matched_good<-read.table(file = "/PATH/TO/fl_meta_matched_good_t.table", header = TRUE)
shared<-fl_meta_matched_good[fl_meta_matched_good$sample_id%in%all$Tumour,]
non_shared<-fl_meta_matched_good[!(fl_meta_matched_good$sample_id%in%all$Tumour),]


##############################

### load cluster summary files for all data sets
gam <- list.files("analysis/fl-dlbcl/mol_time/cluster_summary_gambl_genome--grch37",pattern="*_cluster_summary.txt", full.names = TRUE)

cluster_attack_gam <- lapply(gam,function(x) {
  read.csv(x,  header=TRUE, sep = "\t")
})


### add a column as sample name
for (i in 1:length(cluster_attack_gam)){
  cluster_attack_gam[[i]]<-cbind(cluster_attack_gam[[i]],gam[i])
}

gam1 <- do.call("rbind", cluster_attack_gam) 

####
icgc <- list.files("analysis/fl-dlbcl/mol_time/cluster_summary_icgc_genome--grch37",pattern="*_cluster_summary.txt", full.names = TRUE)

cluster_attack_icgc <- lapply(icgc,function(x) {
  read.csv(x,  header=TRUE, sep = "\t")
})

for (i in 1:length(cluster_attack_icgc)){
  cluster_attack_icgc[[i]]<-cbind(cluster_attack_icgc[[i]],icgc[i])
}


icgc1 <- do.call("rbind", cluster_attack_icgc) 

icgc1$sample<-gsub("analysis/fl-dlbcl/mol_time/cluster_summary_icgc_genome--grch37/", "",gsub("_cluster_summary.txt","",icgc1[,15]))
gam1$sample<-gsub("analysis/fl-dlbcl/mol_time/cluster_summary_gambl_genome--grch37/", "",gsub("_cluster_summary.txt","",gam1[,15]))


icgc1_f<-icgc1[,-15]
gam1_f<-gam1[,-15]


both<-rbind(icgc1_f,gam1_f)
tumour<-data.frame(str_split_fixed(both$sample, "--", 2))
colnames(tumour)<-c("tumour","normal")

both<-cbind(both,tumour) ### mol_time estimates for both ICGC and GAMBL

fl_meta_matched_good<-read.delim(file = "/home/ssoudi/fl_meta_matched_good_all_columns.table", header = TRUE)[,c("patient_id","sample_id","cohort","fl_grade","pathology","consensus_pathology","analysis_cohort","sex")]

#fl_meta_matched_good_t<-fl_meta_matched_good[,c(3,4,5,8,69,70)]
fl_meta_matched_good_t<-fl_meta_matched_good
fl_meta_in_fl_times<-fl_meta_matched_good[fl_meta_matched_good$sample_id%in%both$tumour,]

mol_meta<-merge(both,fl_meta_in_fl_times, by.x = "tumour", by.y ="sample_id")  ### FL+DLBCL

mol_meta$unique_id<-paste0(mol_meta$tumour, "_", mol_meta$code, "_", mol_meta$segment, "_", "_", mol_meta$chrom)


### limit to fl or dlbcl seperately
meta_DLBCL<-fl_meta_matched_good_t[fl_meta_matched_good_t$pathology=="DLBCL",]
meta_FL<-fl_meta_matched_good_t[fl_meta_matched_good_t$pathology=="FL",]

#overlap_with_molID_set_dlbcl<-both[both$tumour%in%meta_DLBCL$sample_id,]
#overlap_with_molID_set_fl<-both[both$tumour%in%meta_FL$sample_id,]

#mol_meta_dlbcl<-merge(overlap_with_molID_set_dlbcl,fl_meta_in_fl_times, by.x = "tumour", by.y ="sample_id")  ### FL+DLBCL
#mol_meta_fl<-merge(overlap_with_molID_set_fl,fl_meta_in_fl_times, by.x = "tumour", by.y ="sample_id")  ### FL+DLBCL

######################################################################
aa<-mol_meta[,c("tumour","segment","plot_dot1")]

colnames(aa)<-c("Sample","Gene","BT_VAF")
#aa$chr <- gsub("-.*$", "", aa$Gene)

qq<-data.frame(str_split_fixed(aa$Gene, "-",4))
#qq[,c(2:3)]<-sapply(qq[,c(2:3)], as.numeric)
colnames(qq)<-c("chr","start","end","seg")

cboth<-cbind(aa,qq)

### read p and q arms start and end coordinate ###
arm<-read.table(file = "/projects/rmorin/projects/gambl-repos/gambl-ssoudi/analysis/mol_time/icgc_dart/icgc_dart/MM_icgc_genome--hs37d5/chromArm.hg19.tsv", header = TRUE)
colnames(arm)<-c("chr","start","end","arm")

chr<-as.character(unique(cboth$chr))



out_res<-NULL
for (i in 1:length(chr)){
  
  cboth_foc<-cboth[cboth$chr==chr[i],]
  arm_foc<-arm[arm$chr==chr[i] &arm$arm == "p" ,]
  arm_foc_q<-arm[arm$chr==chr[i] &arm$arm == "q" ,]
  diff<-abs (arm_foc$end-arm_foc_q$start)
  
  for (j in 1:nrow(cboth_foc)){
    
    cboth_foc_f<-cboth_foc[j,]
    
    if (as.numeric(as.character(cboth_foc_f$start)) > (arm_foc$start) & as.numeric(as.character(cboth_foc_f$end))< (arm_foc$end)+diff-1) {
      cboth_foc_f$id<-paste(cboth_foc_f$chr,arm_foc$arm, sep = "") 
      cboth_foc_f$id<-paste(cboth_foc_f$id,cboth_foc_f$seg, sep = "-") 
      out_res <- rbind (cboth_foc_f,out_res)} else if (as.numeric(as.character(cboth_foc_f$start)) > ((arm_foc_q$start)-diff) & as.numeric(as.character(cboth_foc_f$end))< (arm_foc_q$end)){
        
        
        cboth_foc_f$id<-paste(cboth_foc_f$chr,arm_foc_q$arm, sep = "") 
        cboth_foc_f$id<-paste(cboth_foc_f$id,cboth_foc_f$seg, sep = "-")  
        out_res <- rbind (cboth_foc_f,out_res)} else {
          
          cboth_foc_f$id<-paste(cboth_foc_f$chr,cboth_foc_f$seg, sep = "-") 
          out_res <- rbind (cboth_foc_f,out_res)
        }
    
    #if (as.numeric(as.character(cboth_foc_f$start)) > ((arm_foc_q$start)-diff) & as.numeric(as.character(cboth_foc_f$end))< (arm_foc_q$end)){
    # cboth_foc_f$id<-paste(cboth_foc_f$chr,arm_foc_q$arm, sep = "") 
    #cboth_foc_f$id<-paste(cboth_foc_f$id,cboth_foc_f$seg, sep = "-")  
    #out_res <- rbind (cboth_foc_f,out_res)}
    #if (as.numeric(as.character(cboth_foc_f$start)) < ((arm_foc$end)+diff) & as.numeric(as.character(cboth_foc_f$end))> (arm_foc_q$start)){
    
    #cboth_foc_f$id<-paste(cboth_foc_f$chr,cboth_foc_f$seg, sep = "-") 
    #out_res <- rbind (cboth_foc_f,out_res)
    # }
    
  } ## j loop
  
} ## i loop 


out_res$Gene<-gsub("gain_2","gain",out_res$Gene)
out_res$seg<-gsub("gain_2","gain",out_res$seg)
out_res$id<-gsub("gain_2","gain",out_res$id)

out_res_bed<-out_res[,c("chr" ,"start","end","Sample","BT_VAF","seg","id")]
write.table(out_res_bed, file = "/home/ssoudi/mol_time_intervals.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


################################################################
## laod GISTIC peaks
overlap<-read.table(file = "/home/ssoudi/mol_time_gistic_peaks.table", header = FALSE)       # ovelapped segments with gistic peaks

mol<-data.frame(str_split_fixed(overlap$V11, "-", 2))
mol$X2<-gsub("LOH","Del",mol$X2)
mol$X2<-gsub("gain","Amp",mol$X2)
colnames(mol)<-c("mol_chrom","mol_cnv")
gistic<-data.frame(str_split_fixed(overlap$V4, "_", 2))
colnames(gistic)<-c("gis_chrom","gis_cnv")

overlap_with_molID<-cbind(overlap,mol,gistic)
overlap_with_molID_set<-overlap_with_molID[overlap_with_molID$mol_cnv==overlap_with_molID$gis_cnv,]

#head(overlap_with_molID_set)
#V1      V2      V3          V4   V5      V6       V7        V8        V9 V10       V11 mol_chrom mol_cnv gis_chrom gis_cnv
#1 chr1 3855794 8861261 1p36.23_Del chr1  762601 35427957   SP59300 0.3000000 LOH chr1p-LOH     chr1p     Del   1p36.23     Del
#3 chr1 3855794 8861261 1p36.23_Del chr1  776546 43642673  FL3014T1 0.9620253 LOH chr1p-LOH     chr1p     Del   1p36.23     Del
#4 chr1 3855794 8861261 1p36.23_Del chr1  768253 32597889  FL3003T1 0.4912281 LOH chr1p-LOH     chr1p     Del   1p36.23     Del



### make a dataframe of GISTIC segments and overlapping CNVs
lymphoma_genes<-read.table(file = "/home/ssoudi/lymphoma_genes_within_fl_dlbcl_GISTIC_peaks.table", header = TRUE)

segs<-unique(lymphoma_genes$CNV_peak_id)

out_res<-NULL
for (k in 1:length(segs)){
  
  focal_seg<-lymphoma_genes[lymphoma_genes$CNV_peak_id==segs[k],]
  
  
  if (nrow(focal_seg) ==1){
    
    rr<-focal_seg[,c("CNV_peak_id","lymphoma_gene_id")]
    rr$new_CNV_peak_id<-paste(rr$CNV_peak_id,rr$lymphoma_gene_id, sep = "--")
  }
  
  if (nrow(focal_seg) ==2){
    lymphoma_gene_id<-paste(focal_seg[1,8],focal_seg[2,8],sep = ",")
    lymphoma_gene_id<-as.data.frame(lymphoma_gene_id)
    rr<-cbind(data.frame(focal_seg[1,4]),lymphoma_gene_id)
    colnames(rr)<-c("CNV_peak_id","lymphoma_gene_id")
    rr$new_CNV_peak_id<-paste(rr$CNV_peak_id,rr$lymphoma_gene_id, sep = "--")
    
  }
  
  if (nrow(focal_seg) ==3){
    lymphoma_gene_id<-paste(focal_seg[1,8],focal_seg[2,8],focal_seg[3,8])
    lymphoma_gene_id<-as.data.frame(lymphoma_gene_id)
    rr<-cbind(data.frame(focal_seg[1,4]),lymphoma_gene_id)
    colnames(rr)<-c("CNV_peak_id","lymphoma_gene_id")
    rr$new_CNV_peak_id<-paste(rr$CNV_peak_id,rr$lymphoma_gene_id, sep = "--")
    
  }
  
  out_res<-rbind(rr,out_res)
  #out_res$seg<-paste(out_res$CNV_peak_id, out_res$lymphoma_gene_id, sep = "--")
}

################################
################################


out_res_o<-NULL
for(oo in 1:nrow(overlap_with_molID_set)){
  foc<-overlap_with_molID_set[oo,]
  if (foc$V4 %in%out_res$CNV_peak_id){
    
    foc_out<-unique(out_res[out_res$CNV_peak_id%in%foc$V4,3])
    foc$new_gis_chrom<-foc_out
    
  } else {
    
    foc$new_gis_chrom<-foc$V4
  }
  out_res_o<-rbind(foc,out_res_o)
}


#out_res_dlcbl<-out_res_o[out_res_o$V8%in%mol_meta_dlbcl$tumour,]
#out_res_fl<-out_res_o[out_res_o$V8%in%mol_meta_fl$tumour,]

#######
elli2<-out_res_o[,c(8,16,9)]
#       Sample       Gene    VAF_adj
#1341  SP59448 chrXq-gain 0.08622754
#1340  SP59436 chrXq-gain 0.20689655
#1296  SP59304  chrXq-LOH 0.40157480


colnames(elli2)<- c("Sample","Gene","BT_VAF") ### keep the columns as the original data frame to avoid issues
colnames(elli2)[1:3]<-c("Sample","Gene","BT_VAF")
colnames(elli2)[3]<-"VAF_adj"
elli2$Sample<-as.character(as.character(elli2$Sample))
elli2$Gene<-as.character(as.character(elli2$Gene))

mutationData<- elli2[order(elli2$Sample),]
mutationTable <- (table(elli2[c("Sample","Gene")]) > 0)+0



#mutationData$DP.Tum <- as.numeric(as.character(mutationData$DP.Tum))
vafCn <- as.numeric(as.character(mutationData$VAF_adj))
#vaf <- mutationData$PM.Tum 
mutationData$Sample <- factor(as.character(mutationData$Sample))
vafCn <- as.numeric(as.character(mutationData$VAF_adj))

list_seg<- unique(mutationData$Gene)
clinicalData<-as.data.frame.matrix(table(elli2$Sample,elli2$Gene))

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
        if(vafCn[i] < vafCn[j]){ 
          #Pidgeonhole
          precedence[as.character(mutationData$Gene[i]),as.character(mutationData$Gene[j])] <- precedence[as.character(mutationData$Gene[i]),as.character(mutationData$Gene[j])] + 1
          l <- c(l, list(c(as.character(mutationData$Gene[i]),as.character(mutationData$Gene[j]))))       
        }
      }
    }
  plist[[s]] <- l
}


sample_list<- unique(elli2$sample)

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


# pdf("bradley_terry_integrated_all_mol_time_order_new_copy-number.pdf", width=8, heigh=8)
par(mfrow=c(1,1))
par(mar=c(3,3,3,0),xpd=T )
nCasesGene <- table(factor(unlist(sapply(plist, function(x) unique(unlist(x)))), levels=colnames(precedence)))
w<- which(nCasesGene > 1)  ##### threshold of homemany cases
#w<- which(nCasesGene > 5)  ##### threshold of homemany cases


fit <- btModel(precedence[w,w]+.01) ## Warning: non-integer counts in a binomial glm!
c <- c(0,coef(fit))
names(c) <- colnames(precedence)[w]
o <- rank(c)
v <- pmin(2,sqrt(c(0,diag(chol2inv(fit$qr$qr)))))
l <- names(c)
m <- paste("n=",nCasesGene[w], sep="")

### for coloring
wnam<-names(nCasesGene[w])
wnam1<-str_split_fixed(wnam, "--",2)[,1]

str_split_fixed(wnam1, "_",2)[,2]

pdf("/home/ssoudi/FL--pathology_MolTime_aggregated_gistic_peaks_annotatted.pdf", width=8, heigh=8)
#options(bitmapType='cairo')
#png("/home/ssoudi/fl_DLBCL_MolTime_aggregated_gistic_peaks_annotatted.png", width=350, heigh=350)
#jpeg("/home/ssoudi/fl_DLBCL_MolTime_aggregated_gistic_peaks_annotatted.jpg", width = 800, height = 800)


par(mfrow=c(1,1))
par(mar=c(4,6,3,0),xpd=T )
plot(-c, o, xlab="Relative time", yaxt="n", pch=19,  ylab="", xlim=c(-5,15), bty="n", main = "FL")
segments(-c-v, o,-c+v,o, col="grey")
par(new=TRUE)
plot(-c, o, xlab="", yaxt="n",xaxt="n", pch=19, 
     ylab="", xlim=c(-5,15), bty="n")
text(-c-v ,o,gsub("_"," ",l), font=3, pos=2, cex=0.8,col = ifelse(str_split_fixed(wnam1, "_",2)[,2]=="Amp","forestgreen","blueviolet"))
text(-c+v ,o,m, font=1, pos=4, cex=1)
dev.off()



