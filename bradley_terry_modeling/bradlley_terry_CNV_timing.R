### This script takes the results of molecular time (mol_time or muttaiontimer) to estimate the order of aquisition of CNVs based omn Bradley-Terry Model.
### To learn about the aproach see this article: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1476-3

## author: Shaghayegh Soudi
## shaghayegh.soudi@gmail.com
# rm(list = ls())


library(stringr)


### list and load input files 

## load cluster_summary files
filenames <- list.files("cluster_summary_multiple_myeloma",pattern="*.txt", full.names = TRUE)
attackStats <- lapply(filenames,function(x) {
  read.csv(x,  header=TRUE, sep = "\t")[,c(1:3)]
})

### add sample name as a column
for (i in 1:length(attackStats)){
  attackStats[[i]]<-cbind(attackStats[[i]],filenames[i])
}
aa <- do.call("rbind", attackStats) 
aa$Sample<-gsub("cluster_summary_multiple_myeloma/", "",gsub("_cluster_summary.txt","",aa[,4]))
aa<-aa[,c(5,2,3)

colnames(aa)<-c("Sample","Gene","BT_VAF")
#aa$chr <- gsub("-.*$", "", aa$Gene)

qq<-data.frame(str_split_fixed(aa$Gene, "-",4))
#qq[,c(2:3)]<-sapply(qq[,c(2:3)], as.numeric)
colnames(qq)<-c("chr","start","end","seg")

cboth<-cbind(aa,qq)

#####
## laod chromosomal coordinate information based on the genome assembly chosen (in this example data mappled on hg19)
arm<-read.table(file = "PATH/TO/chromArm.hg19.tsv", header = TRUE)
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
  
#######
elli2<-out_res[,c(1,8,3)]

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
w<- which(nCasesGene > 1)
fit <- btModel(precedence[w,w]+.01) ## Warning: non-integer counts in a binomial glm!
c <- c(0,coef(fit))
names(c) <- colnames(precedence)[w]
o <- rank(c)
v <- pmin(2,sqrt(c(0,diag(chol2inv(fit$qr$qr)))))
l <- names(c)
m <- paste("n=",nCasesGene[w], sep="")

### for coloring
wnam<-names(nCasesGene[w])
str_split_fixed(wnam, "-",2)[,2]

pdf("cluster_summary/bradley_terry_ICGC_multiplemyeloma_threshold_50.pdf", width=8, heigh=10)
par(mfrow=c(1,1))
#par(mfrow=c(1,2), mai = c(0.1, 0.1, 0.1, 6))
par(mar=c(4,6,3,0),xpd=T )
#par(mar=c(4,3,3,0),xpd=T )
plot(-c, o, xlab="Relative time", yaxt="n", pch=19,  ylab="", xlim=c(-10,10), bty="n")
segments(-c-v, o,-c+v,o, col="grey")
par(new=TRUE)
plot(-c, o, xlab="", yaxt="n",xaxt="n", pch=19, 
     ylab="", xlim=c(-10,10), bty="n", main = "ordering diagram for FLs-GenomeCanada")
text(-c-v ,o,gsub("_"," ",l), font=3, pos=2, cex=0.8,col = ifelse(str_split_fixed(wnam, "-",2)[,2]=="gain","forestgreen","blueviolet"))
text(-c+v ,o,m, font=1, pos=4, cex=1)
dev.off()


