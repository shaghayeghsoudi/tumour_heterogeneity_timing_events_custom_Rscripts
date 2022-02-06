#!/usr/bin/env Rscript

## rm(list = ls())
## this script times CNV segments from molt_time or muttaiontimeR packages and looks at their overlaps with GISTIC peaks and incorporate the involved genes in the timing plot


##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 

library(stringr)

## merge all results from different genome assemblies if applicable
grch37<-read.table(file = "PATH/TO/sample_normal_contamination.txt", header = TRUE)
grch37_TN<-as.data.frame(str_split_fixed(grch37$SAMPLE, "--", 2))
colnames(grch37_TN)<-c("Tumour","Normal")
all<-cbind(grch37,grch37_TN)

#all<-rbind(icgc_grch37,gambl_grch37)

## metadata file for your study 
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


####### make pie charts ########

#example: head(timedGains_PCAWG)
#samplename	chr	start	end	length	time	time_lCI	time_uCI	no.snvs	type	time_lCI_adj	time_uCI_adj	histology_abbreviation
#0009b464-b376-4fbc-8a56-da538269a02f	WGD	NA	NA	NA	0.254159807	NA	NA	6364	WGD	NA	NA	Ovary-AdenoCA
#0040b1b6-b07a-4b6e-90ef-133523eaf412	1	106135000	121500000	15365077	0.001598973	3.50E-05	0.192132294	76	SingleGain	3.28E-05	0.242000671	Liver-HCC


#### deletions #####
summaryTable = fl_meta_matched_good[,c("sample_id","pathology")]
colnames(summaryTable)<-c("samplename","histology_abbreviation")


timedGains = overlap_with_molID_set[,c(1,2,3,4,6,7,8,9,10,11,12,13,14,15)]
timedGains<-timedGains[timedGains$mol_cnv=="Del",]

colnames(timedGains)<-c("chr","start","end","seg_ID","start_seg","end_seg","samplename","time","original_seg","mol_cnv","gis_chrom","type","seg","cnv")

mol_meta_good<-fl_meta_matched_good[,c("sample_id","pathology")]
timedGains<-merge(timedGains,mol_meta_good, by.x = "samplename", by.y ="sample_id")
timedGains$chr<-gsub("chr","",timedGains$chr)

# Add segment identifier
timedGains$segId = paste0(timedGains$samplename, "_", timedGains$chr, "_", timedGains$start, "_", timedGains$end, "_", timedGains$type)
colnames(timedGains)[15]<-"histology_abbreviation"

# Add an "all" chromosome set of rows
#timedGains.all = timedGains
#timedGains.all$chr = "all"
#timedGains = rbind(timedGains, timedGains.all)

# Add n for each cancer type, all results
n.samples.hist = table(unique(timedGains[,which(colnames(timedGains) %in% c("samplename", "histology_abbreviation"))])$histology_abbreviation)
timedGains$n_samples_hist = n.samples.hist[match(timedGains$histology_abbreviation, names(n.samples.hist))]
timedGains$n_histology_abbreviation = paste0(timedGains$histology_abbreviation, " (n = ", timedGains$n_samples_hist, ")")

tumour_types<-sort(unique(timedGains$histology_abbreviation))
tumour_types<-c("DLBCL","FL")

seg_type<-unique(timedGains$seg_ID)
#[1] "10q26.3_Del"                    "11q23.2_Amp--ATM,POU2AF1"      
#[3] "12p13.2_Del"                    "12q12_Amp"                     


#l_del<-l[grep("Del",l)]
#l_del<-gsub("--.*", "", l_del)
#l_del<-gsub("_Del","",l_del)

l_del<-seg_type[grep("Del",seg_type)]
l_del<-gsub("--.*", "", l_del)
l_del<-gsub("_Del","",l_del)


########
# Get proportion of events, per chromosome per cancer type
for (c in 1:length(l_del)){
  #for (h in sort(unique(timedGains$histology_abbreviation))){
  
  for (h in 1:length(tumour_types)) {
    timedGains_CH = timedGains[which(timedGains$seg==l_del[c] & timedGains$histology_abbreviation==tumour_types[h]),]
    #timedGains_CH = timedGains_CH[!(duplicated(timedGains_CH$samplename)),]
    
    
    # Define proportion of events - this is the size of the pie chart
    n_samples = length(unique(timedGains_CH$sample))
    n_samples_total = length(subset(summaryTable, histology_abbreviation==tumour_types[h])$samplename)
    prop = sqrt((n_samples/n_samples_total)/pi)
    
    #timedGains<-timedGains[timedGains$seg_ID==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h],]
    timedGains$prop[timedGains$seg==l_del[c] & timedGains$histology_abbreviation==tumour_types[h]] = prop
    
    # Define start/end boundaries for segments
    timedGains_CH = timedGains_CH[order(-timedGains_CH$time),]
    timedGains_CH$start_segment = cumsum(rep(1/nrow(timedGains_CH), nrow(timedGains_CH))) - 1/nrow(timedGains_CH)
    timedGains_CH$end_segment = cumsum(rep(1/nrow(timedGains_CH), nrow(timedGains_CH)))
    
    timedGains$start_segment[timedGains$seg==l_del[c] & timedGains$histology_abbreviation==tumour_types[h]] = timedGains_CH$start_segment[match(timedGains$segId[timedGains$seg==l_del[c] & timedGains$histology_abbreviation==tumour_types[h]], timedGains_CH$segId)]
    
    timedGains$end_segment[timedGains$seg==l_del[c] & timedGains$histology_abbreviation==tumour_types[h]] = timedGains_CH$end_segment[match(timedGains$segId[timedGains$seg==l_del[c] & timedGains$histology_abbreviation==tumour_types[h]], timedGains_CH$segId)]
    
  }
}


timedGains<-timedGains[!is.na(timedGains$prop),]



# Sort cancer types by median pi0
medianTime = aggregate(timedGains$time, list(timedGains$n_histology_abbreviation), median)
timedGains$n_histology_abbreviation = factor(timedGains$n_histology_abbreviation, medianTime[order(medianTime$x),]$Group.1)
# Remove cancer types with fewer than 100 segments, 15 samples, and mean pie size < 0.2
#timedGains = subset(timedGains, histology_abbreviation %in% names(which(table(timedGains$histology_abbreviation) > 100)))
quiet.ttypes = c(names(which(table(summaryTable$histology_abbreviation) < 20)))
mean_sizes = aggregate(timedGains$prop[timedGains$chr!="all"], list(timedGains$histology_abbreviation[timedGains$chr!="all"]), mean)
quiet.ttypes = unique(c(quiet.ttypes,  as.character(mean_sizes[which(mean_sizes$x < 0.12),]$Group.1)))



Plot as pie charts
```{r fig1, fig.height=22, fig.width=25}

theme_set(theme_grey())
ggplot(data=timedGains[!timedGains$histology_abbreviation %in% quiet.ttypes,]) +
  #facet_grid(n_histology_abbreviation~factor(chr, levels=c(1:22, "X", "all", "WGD"))) + 
  facet_grid(n_histology_abbreviation~factor(seg, levels=l_del)) + 
  
  
  
  coord_polar(theta="x") + 
  geom_rect(aes(xmin=start_segment, xmax=end_segment, ymin=0, ymax=prop, fill=time)) +
  #geom_rect(aes(xmin=0, xmax=0.1, ymin=0, ymax=prop, fill=time)) +
  
  
  #scale_fill_distiller(palette = "PRGn", name="Molecular Time", breaks=c(0,0.25,0.5,0.75,1,1.25,1.5), limits=c(0,1.5)) +
  scale_fill_distiller(palette = "RdBu", name="Molecular Time", breaks=c(0,0.5,1,1.5), limits=c(0,1.5)) +
  
  theme(strip.background = element_blank()) +
  theme(strip.text.y=element_text(angle=0, size=15, hjust=0)) + 
  theme(strip.text.x=element_text(angle=17,size=11)) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(), panel.grid=element_blank()) +      ## panel 
  theme(panel.spacing=unit(0, "lines"),
        panel.border=element_rect(color="gray", fill=NA, size=0.8)) +
  theme(panel.background=element_blank()) +  ###
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) + 
  theme(legend.position="right") +
  guides(fill = guide_colorbar(barwidth=1.5, barheight=8)) +
  ggtitle("Deletions mutational timing estimates") +
  theme(plot.title = element_text(hjust = 0.5))

```



##############################
###### Amplifications ########
summaryTable = fl_meta_matched_good[,c("sample_id","consensus_pathology")]
colnames(summaryTable)<-c("samplename","histology_abbreviation")


timedGains = overlap_with_molID_set[,c(1,2,3,4,6,7,8,9,10,11,12,13,14,15)]
timedGains<-timedGains[timedGains$mol_cnv=="Amp",]

colnames(timedGains)<-c("chr","start","end","seg_ID","start_seg","end_seg","samplename","time","original_seg","mol_cnv","gis_chrom","type","seg","cnv")

mol_meta_good<-fl_meta_matched_good[,c("sample_id","consensus_pathology")]
timedGains<-merge(timedGains,mol_meta_good, by.x = "samplename", by.y ="sample_id")
timedGains$chr<-gsub("chr","",timedGains$chr)

# Add segment identifier
timedGains$segId = paste0(timedGains$samplename, "_", timedGains$chr, "_", timedGains$start, "_", timedGains$end, "_", timedGains$type)
colnames(timedGains)[15]<-"histology_abbreviation"

# Add an "all" chromosome set of rows
#timedGains.all = timedGains
#timedGains.all$chr = "all"
#timedGains = rbind(timedGains, timedGains.all)

# Add n for each cancer type, all results
n.samples.hist = table(unique(timedGains[,which(colnames(timedGains) %in% c("samplename", "histology_abbreviation"))])$histology_abbreviation)
timedGains$n_samples_hist = n.samples.hist[match(timedGains$histology_abbreviation, names(n.samples.hist))]
timedGains$n_histology_abbreviation = paste0(timedGains$histology_abbreviation, " (n = ", timedGains$n_samples_hist, ")")

tumour_types<-sort(unique(timedGains$histology_abbreviation))
tumour_types<-c("DLBCL","FL")

seg_type<-unique(timedGains$seg_ID)

#l
#[1] "10q26.3_Del"                    "11q23.2_Amp--ATM,POU2AF1"      
#[3] "12p13.2_Del"                    "12q12_Amp"                     
#[5] "12q15_Amp"                      "13q12.3_Del" 

l_amp<-seg_type[grep("Amp",seg_type)]
l_amp<-gsub("--.*", "", l_amp)
l_amp<-gsub("_Amp","",l_amp)


##
# Get proportion of events, per chromosome per cancer type
for (c in 1:length(l_amp)){
  #for (h in sort(unique(timedGains$histology_abbreviation))){
  
  for (h in 1:length(tumour_types)) {
    timedGains_CH = timedGains[which(timedGains$seg==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h]),]
    
    # Define proportion of events - this is the size of the pie chart
    n_samples = length(unique(timedGains_CH$sample))
    n_samples_total = length(subset(summaryTable, histology_abbreviation==tumour_types[h])$samplename)
    prop = sqrt((n_samples/n_samples_total)/pi)
    
    #timedGains<-timedGains[timedGains$seg_ID==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h],]
    timedGains$prop[timedGains$seg==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h]] = prop
    
    # Define start/end boundaries for segments
    timedGains_CH = timedGains_CH[order(-timedGains_CH$time),]
    timedGains_CH$start_segment = cumsum(rep(1/nrow(timedGains_CH), nrow(timedGains_CH))) - 1/nrow(timedGains_CH)
    timedGains_CH$end_segment = cumsum(rep(1/nrow(timedGains_CH), nrow(timedGains_CH)))
    
    timedGains$start_segment[timedGains$seg==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h]] = timedGains_CH$start_segment[match(timedGains$segId[timedGains$seg==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h]], timedGains_CH$segId)]
    
    timedGains$end_segment[timedGains$seg==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h]] = timedGains_CH$end_segment[match(timedGains$segId[timedGains$seg==l_amp[c] & timedGains$histology_abbreviation==tumour_types[h]], timedGains_CH$segId)]
    
  }
}


timedGains<-timedGains[!is.na(timedGains$prop),]



# Sort cancer types by median pi0
medianTime = aggregate(timedGains$time, list(timedGains$n_histology_abbreviation), median)
timedGains$n_histology_abbreviation = factor(timedGains$n_histology_abbreviation, medianTime[order(medianTime$x),]$Group.1)
# Remove cancer types with fewer than 100 segments, 15 samples, and mean pie size < 0.2
#timedGains = subset(timedGains, histology_abbreviation %in% names(which(table(timedGains$histology_abbreviation) > 100)))
quiet.ttypes = c(names(which(table(summaryTable$histology_abbreviation) < 20)))
mean_sizes = aggregate(timedGains$prop[timedGains$chr!="all"], list(timedGains$histology_abbreviation[timedGains$chr!="all"]), mean)
quiet.ttypes = unique(c(quiet.ttypes,  as.character(mean_sizes[which(mean_sizes$x < 0.12),]$Group.1)))
```


Plot as pie charts
```{r fig1, fig.height=22, fig.width=25}

theme_set(theme_grey())
ggplot(data=timedGains[!timedGains$histology_abbreviation %in% quiet.ttypes,]) +
  #facet_grid(n_histology_abbreviation~factor(chr, levels=c(1:22, "X", "all", "WGD"))) + 
  facet_grid(n_histology_abbreviation~factor(seg, levels=l_amp)) + 
  
  
  
  coord_polar(theta="x") + 
  geom_rect(aes(xmin=start_segment, xmax=end_segment, ymin=0, ymax=prop, fill=time)) +
  #geom_rect(aes(xmin=0, xmax=0.1, ymin=0, ymax=prop, fill=time)) +
  
  
  scale_fill_distiller(palette = "PRGn", name="Molecular Time", breaks=c(0,0.25,0.5,0.75,1,1.25,1.5), limits=c(0,1.5)) +
  #scale_fill_distiller(palette = "RdBu", n=6, name="Time, fraction\n of mutations", breaks=c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6), limits=c(0,1.6)) +
  
  theme(strip.background = element_blank()) +
  theme(strip.text.y=element_text(angle=0, size=15, hjust=0)) + 
  theme(strip.text.x=element_text(angle=17,size=14)) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(), panel.grid=element_blank()) +      ## panel 
  theme(panel.spacing=unit(0, "lines"),
        panel.border=element_rect(color="gray", fill=NA, size=0.8)) +
  theme(panel.background=element_blank()) +  ###
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) + 
  theme(legend.position="right") +
  guides(fill = guide_colorbar(barwidth=1.5, barheight=8)) +
  ggtitle("Amplifications mutational timing estimates") +
  theme(plot.title = element_text(hjust = 0.5))

```


###############################################
###### complex-heatmaps (moleculat time) ######
###############################################
#rm(list = ls())
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

setwd("/projects/rmorin/projects/gambl-repos/gambl-ssoudi/")


fl_meta_matched_good<-read.delim(file = "/home/ssoudi/fl_meta_matched_good_all_columns.table", header = TRUE)[,c("patient_id","sample_id","cohort","fl_grade","pathology","consensus_pathology","analysis_cohort","sex")]

#fl_meta_matched_good_t<-fl_meta_matched_good[,c(3,4,5,8,69,70)]
fl_meta_matched_good_t<-fl_meta_matched_good



### load gambl cluster_summery files
gam <- list.files("analysis/fl-dlbcl/mol_time/cluster_summary_gambl_genome--grch37",pattern="*_cluster_summary.txt", full.names = TRUE)

cluster_attack_gam <- lapply(gam,function(x) {
  read.csv(x,  header=TRUE, sep = "\t")
})


### add a column as sample name
for (i in 1:length(cluster_attack_gam)){
  cluster_attack_gam[[i]]<-cbind(cluster_attack_gam[[i]],gam[i])
}
gam1 <- do.call("rbind", cluster_attack_gam) 

### load icgc cluster_summery files
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

## merge gambl and icgc together
both<-rbind(icgc1_f,gam1_f)
tumour<-data.frame(str_split_fixed(both$sample, "--", 2))
colnames(tumour)<-c("tumour","normal")
both<-cbind(both,tumour)

#fl_meta_matched_good_t<-fl_meta_matched_good[,c(3,4,5,8,69,70)]
fl_meta_matched_good_t<-fl_meta_matched_good
fl_meta_in_fl_times<-fl_meta_matched_good[fl_meta_matched_good$sample_id%in%both$tumour,]
mol_meta<-merge(both,fl_meta_matched_good_t[,c("sample_id","cohort","pathology","analysis_cohort","sex")], by.x = "tumour", by.y = "sample_id")


overlap<-read.table(file = "/home/ssoudi/mol_time_gistic_peaks.table", header = FALSE)       # ovelapped segments with gistic peaks

mol<-data.frame(str_split_fixed(overlap$V11, "-", 2))
mol$X2<-gsub("LOH","Del",mol$X2)
mol$X2<-gsub("gain","Amp",mol$X2)
colnames(mol)<-c("mol_chrom","mol_cnv")
gistic<-data.frame(str_split_fixed(overlap$V4, "_", 2))
colnames(gistic)<-c("gis_chrom","gis_cnv")

overlap_with_molID<-cbind(overlap,mol,gistic)
overlap_with_molID_set<-overlap_with_molID[overlap_with_molID$mol_cnv==overlap_with_molID$gis_cnv,]
colnames(overlap_with_molID_set)<-c("chrom","g_start","g_end","g_ID","chrom2","sart","end","sample","time","seg","seg_ID","mol_chrom", "mol_cnv", "gis_chrom", "gis_cnv")

#head(overlap_with_molID_set)
#V1      V2      V3          V4   V5      V6       V7        V8        V9 V10       V11 mol_chrom mol_cnv gis_chrom gis_cnv
#1 chr1 3855794 8861261 1p36.23_Del chr1  762601 35427957   SP59300 0.3000000 LOH chr1p-LOH     chr1p     Del   1p36.23     Del
#3 chr1 3855794 8861261 1p36.23_Del chr1  776546 43642673  FL3014T1 0.9620253 LOH chr1p-LOH     chr1p     Del   1p36.23     Del
#4 chr1 3855794 8861261 1p36.23_Del chr1  768253 32597889  FL3003T1 0.4912281 LOH chr1p-LOH     chr1p     Del   1p36.23     Del

moltime_ovelpa_gistic_samples<-unique(overlap_with_molID_set$sample)


out_res<-NULL
for (ss in 1:length(moltime_ovelpa_gistic_samples)){
  
  #for (ss in 1:54){ 
  
  foc<-t(overlap_with_molID_set[overlap_with_molID_set$sample==moltime_ovelpa_gistic_samples[ss],c("gis_chrom","time","g_ID")])
  
  cl.names <- foc["g_ID",]
  colnames(foc) <- cl.names
  
  if(ncol(foc)> 1){
    
    foc_r<-as.data.frame(foc)["time",]
    foc_r$sample<-moltime_ovelpa_gistic_samples[ss]
    rownames(foc_r)<-NULL
    
    meta_foc<-fl_meta_matched_good_t[fl_meta_matched_good_t$sample_id==as.character(moltime_ovelpa_gistic_samples[ss]),c("cohort","pathology","fl_grade","analysis_cohort","sex")]
    
    foc_mol_meta<-cbind(foc_r,meta_foc)
    
  }## if1
  
       if (ncol(foc)== 1){
  
         foc_r<-data.frame((foc)["time",])
         colnames(foc_r)<-foc["g_ID",]
         
         foc_r$sample<-moltime_ovelpa_gistic_samples[ss]
         rownames(foc_r)<-NULL
         
         meta_foc<-fl_meta_matched_good_t[fl_meta_matched_good_t$sample_id==as.character(moltime_ovelpa_gistic_samples[ss]),c("cohort","pathology","fl_grade","analysis_cohort","sex")]
         
         foc_mol_meta<-cbind(foc_r,meta_foc)
  
       } ## if2

  out_res<-rbind.fill(foc_mol_meta,out_res)
}

rownames(out_res)<-out_res[,"sample"]
out_res<-out_res[order(out_res$analysis_cohort),]

###
### out_res_del
names<-colnames(out_res)
names_de<-grep("Del",names)

out_res_del<-out_res[,c(2,3,4,5,6,7,names_de)]
out_res_del<-out_res_del[!(rowSums(is.na(out_res_del[,7:length(out_res_del)]))>35),]


out_res_del_mat<-out_res_del[,(grep("Del",names(out_res_del)))]



i <- c(1:ncol(out_res_del_mat))  
out_res_del_mat<-out_res_del_mat[ , i] <- apply(out_res_del_mat[ , i], 2,            # Specify own function within apply
                                        function(x) as.numeric(as.character(x)))


colnames(out_res_del_mat)<-gsub("_Del","",colnames(out_res_del_mat))
col = colorRamp2(c(1.5,1,0.5,0), brewer.pal(n=4, name="RdBu"))

pdf(file = "heatmap_FL-DLBCL_MolTime_gistic_peaks_Delitions.pdf", width=10, heigh=10)
Heatmap(out_res_del_mat, name = "molecular time",  na_col = "gray97",
        show_row_names = FALSE, rect_gp = gpar(col = "white", lwd = 2),
        column_title = "FL-DLBCL Delitions Molecular Time",width = ncol(out_res_del_mat)*unit(5, "mm"), 
        #height = nrow(mat)*unit(4, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),column_names_rot = 45,
        show_column_names = TRUE,cluster_columns = FALSE, cluster_rows = FALSE, col = col,heatmap_legend_param = list(legend_height = unit(3, "cm")))+
  Heatmap(out_res_del$analysis_cohort, name = "analysis_cohort", width = unit(5, "mm"), c("COM" = "darkorange3", "DLBCL" = "darkorange","FL" = "coral2", "tFL" = "darkolivegreen4"))+
  Heatmap(out_res_del$pathology, name = "pathology", width = unit(5, "mm"), c("DLBCL" = "firebrick4", "FL" = "firebrick3"))+
  Heatmap(out_res_del$fl_grade, name = "fl_grade", width = unit(5, "mm"), c("NA" = "gray60", "FOLL1" = "darkorange","FOLL3A" = "coral2", "FOLL2" = "darkolivegreen4"))+
  Heatmap(out_res_del$sex, name = "sex", width = unit(5, "mm"), c("NA" = "gray60", "M" = "steelblue4","F" = "tomato2"))
dev.off()  

##
## deal with NA

#ee<-out_res_del_mat[rowSums(is.na(out_res_del_mat)) == ncol(out_res_del_mat), ]
#rr<-out_res_del_mat[rowSums(is.na(out_res_del_mat)) < 30,]

Heatmap(rr)










#################
### out_res_amp
names<-colnames(out_res)
names_amp<-grep("Amp",names)

out_res_amp<-out_res[,c(2,3,4,5,6,7,names_amp)]
out_res_amp<-out_res_amp[!(rowSums(is.na(out_res_amp[,7:length(out_res_amp)]))>18),]


out_res_amp_mat<-out_res_amp[,(grep("Amp",names(out_res_amp)))]



i <- c(1:ncol(out_res_amp_mat))  
out_res_amp_mat<-out_res_amp_mat[ , i] <- apply(out_res_amp_mat[ , i], 2,            # Specify own function within apply
                                                function(x) as.numeric(as.character(x)))


colnames(out_res_amp_mat)<-gsub("_Amp","",colnames(out_res_amp_mat))
col = colorRamp2(c(1.5,1,0.5,0), brewer.pal(n=4, name="RdBu"))

pdf(file = "heatmap_FL-DLBCL_MolTime_gistic_peaks_Amplifications.pdf", width=10, heigh=10)
Heatmap(out_res_amp_mat, name = "molecular time",  na_col = "gray97",
        show_row_names = FALSE, rect_gp = gpar(col = "white", lwd = 2),
        column_title = "FL-DLBCL Amplifications Molecular Time",width = ncol(out_res_del_mat)*unit(5, "mm"), 
        #height = nrow(mat)*unit(4, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),column_names_rot = 45,
        show_column_names = TRUE,cluster_columns = FALSE, cluster_rows = FALSE, col = col,heatmap_legend_param = list(legend_height = unit(3, "cm")))+
  Heatmap(out_res_amp$analysis_cohort, name = "analysis_cohort", width = unit(5, "mm"), c("COM" = "darkorange3", "DLBCL" = "darkorange","FL" = "coral2", "tFL" = "darkolivegreen4"))+
  Heatmap(out_res_amp$pathology, name = "pathology", width = unit(5, "mm"), c("DLBCL" = "firebrick4", "FL" = "firebrick3"))+
  Heatmap(out_res_amp$fl_grade, name = "fl_grade", width = unit(5, "mm"), c("NA" = "gray60", "FOLL1" = "darkorange","FOLL3A" = "coral2", "FOLL2" = "darkolivegreen4"))+
  Heatmap(out_res_amp$sex, name = "sex", width = unit(5, "mm"), c("NA" = "gray60", "M" = "steelblue4","F" = "tomato2"))
dev.off()  


