#!/usr/bin/env Rscript

#rm(list = ls())
## this script convert battenberg subclonall file CNV file format (i.e., copy_number.txt) required by mol_time package


##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 


##### convert battenberg subclonal to mol_time ######

### example input from mol_time documenttaion
#sample	Chrom	start	end	class	nMaj1_A	nMin1_A	frac1_A	nMaj2_A	nMin2_A	code_batt
#PD26400a	1	762601	84710840	normal	1	1	1	NA	NA	no_sub
#PD26400a	1	84714419	94909862	deletion	1	0	1	NA	NA	no_sub


### list and load input files 
filenames <- list.files("genome--hg19",pattern="*_subclones.txt", full.names = TRUE)
attackStats <- lapply(filenames,function(x) {
      read.csv(x,  header=TRUE, sep = "\t")[,c(1:4,8:12)]
    })

### add sample name as a column
for (i in 1:length(attackStats)){
    attackStats[[i]]<-cbind(attackStats[[i]],filenames[i])
    }
aa <- do.call("rbind", attackStats) 
#head(aa)

aa$chr<-gsub("chr","",aa$chr)
aa$zz<-gsub("genome--hg19/", "",gsub("_subclones.txt","",aa[,10]))

### make a column to assign clonal and sub_clonal status
aa$class = ifelse(is.na(aa$nMaj2_A) ,"no_sub","subclonal")

### added March 10 ####

### make a temporary column
aa$qq<- ifelse (aa$class== "subclonal",
    ifelse (aa$class== "subclonal" & aa$frac1_A>=0.9 , "subclonal>90",
    ifelse(aa$frac1_A<=0.09 , "subclona<10","subclonal")),"subclonal")
    
aa$class<-ifelse(aa$qq== "subclonal>90", aa$qq, aa$class)
aa$class<-ifelse(aa$qq== "subclona<10", aa$qq, aa$class)


### make a column to assign normal, gain, deletion, LOH, WCD status
aa$code_batt<- ifelse(aa$class== "no_sub" | aa$class=="subclonal>90", 
    ifelse (aa$nMaj1_A==1 & aa$nMin1_A==1, "normal",
    ifelse(aa$nMaj1_A<=1 & aa$nMin1_A==0, "deletion",
    ifelse(aa$nMaj1_A==2 & aa$nMin1_A==1, "gain",
    ifelse(aa$nMaj1_A==3 & aa$nMin1_A==1, "gain_2",
    ifelse(aa$nMaj1_A==4 & aa$nMin1_A==1, "gain_3", 
    ifelse(aa$nMaj1_A>=2 & aa$nMin1_A>=2, "Whole_chromosome_duplication",
    ifelse(aa$nMaj1_A==2 & aa$nMin1_A==0, "LOH",
    ifelse(aa$nMaj1_A==3 & aa$nMin1_A==0, "LOH_2",
    ifelse(aa$nMaj1_A==4 & aa$nMin1_A==0, "LOH_3",
    ifelse(aa$nMaj1_A>=5 & aa$nMin1_A==0, "LOH_High","-")
    ))))))))),"-")
 

aa$fake<- ifelse(aa$qq== "subclona<10",
    ifelse (aa$nMaj2_A==1 & aa$nMin2_A==1, "normal",
    ifelse(aa$nMaj2_A<=1 & aa$nMin2_A==0, "deletion",
    ifelse(aa$nMaj2_A==2 & aa$nMin2_A==1, "gain",
    ifelse(aa$nMaj2_A==2 & aa$nMin2_A==0, "LOH","-")))),"-")


aa$code_batt<-ifelse(aa$qq== "subclona<10", aa$fake, aa$code_batt)    


## adjust colnames and order colmns according to the maunal
aa_fix<-aa[!(aa$nMaj1_A==0 & aa$nMin1_A==0),] 
aa<-aa[,c(11,1,2,3,14,5:9,12)]
colnames(aa)<-c("sample","Chrom","start","end","class","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","code_batt")
write.table(aa, file = "/PATH/TO/YOUR/DIR/copy_number.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")

