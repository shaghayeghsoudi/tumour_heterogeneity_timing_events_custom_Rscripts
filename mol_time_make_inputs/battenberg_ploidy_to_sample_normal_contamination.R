#!/usr/bin/env Rscript

#rm(list = ls())
## this script convert battenberg ploidy file to sample_normal_contamination.txt file format required by mol_time package


##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 


##### sample_normal_contamination.txt ######

#SAMPLE	ABBR_CELL_FRAC	TUM_PLOIDY
#PD26400a	0.94	2.205113821
#PD26400c	1	2.169461186


### list and load input files 
purity <- list.files("genome--hg19/", pattern="*cellularity_ploidy.txt", full.names = TRUE)
attackpurity <- lapply(purity,function(x) {
      read.csv(x,  header=TRUE, sep = "\t")[,c(1:2)]
    })

### add sample name
for (i in 1:length(attackpurity)){
    attackpurity[[i]]<-cbind(attackpurity[[i]],purity[i])
    }
purity_table <- do.call("rbind", attackpurity) 

purity_table$SAMPLE<-gsub("genome--hg19/", "",gsub("_cellularity_ploidy.txt","",purity_table[,3]))

purity_table<-purity_table[,c(4,1,2)]
colnames(purity_table)<-c("SAMPLE","ABBR_CELL_FRAC","TUM_PLOIDY")

write.table(purity_table,file = "/projects/rmorin/projects/gambl-repos/gambl-ssoudi/analysis/mol_time/gambl/genome--hg38/data_genome--grch37/sample_normal_contamination.txt" ,col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")

