library(ggplot2)
library(plyr)


file.names<-list.files("variant_type_counts", pattern  = "*.table", full.names = TRUE)

cluster_attack<-lapply(file.names, function(x){
  read.delim(x,  header=TRUE, sep = "\t")
})


### add a column as sample name

#out_res<-NULL
#for (i in 1:length(cluster_attack)){
#  
#  ww<-cluster_attack[[i]][,c(1,2)]
#  ww_dat<-data.frame(t(ww),row.names = NULL)
#  colnames(ww_dat)<-ww_dat[1,]
#  ww_dat<-ww_dat[-1,]
#  
#  ww_dat$sample<-file.names[i]
#  ww_dat$Sample<-gsub("variant_type_counts/", "",gsub("_variant_type_counts.table","",ww_dat$sample))
#  out_res<-rbind(ww_dat,out_res)
#  
#  }
#

#out_res_final<-out_res[,c("Sample","clonalearly","clonallate","clonalNA","subclonal")]
#barplot(out_res_final,beside = FALSE)



for (i in 1:length(cluster_attack)){
  cluster_attack[[i]]<-cbind(cluster_attack[[i]],file.names[i])
}

aa <- do.call("rbind", cluster_attack) 
aa$Sample<-gsub("variant_type_counts/", "",gsub("_variant_type_counts.table","",aa[,3]))
aa<-aa[,c("Sample","variant_type","frequency")]


ggplot(aa, aes(fill=variant_type, y=frequency, x=Sample)) + 
  geom_bar(position='stack', stat='identity') +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 3.5))


######
##### count driver clonal and subclonal mutations and plot it 
file.names<-list.files("variants_annotatted", pattern  = "*.table", full.names = TRUE)

cluster_attack<-lapply(file.names, function(x){
  read.delim(x,  header=TRUE, sep = "\t")
})



for (i in 1:length(cluster_attack)){
  cluster_attack[[i]]<-cbind(cluster_attack[[i]],file.names[i])
}

aa <- do.call("rbind", cluster_attack) 
aa$Sample<-gsub("variants_annotatted/", "",gsub("_variants_annotatted.table","",aa[,20]))
aa<-aa[,c("Sample","chrom","pos","CLS")]

## load gene coordinate
gene_pos<-read.delim(file = "~/Dropbox/cancer_reserach/Genome_canada_report/cll_genes.txt", header = TRUE)
genes<-unique(gene_pos$name)
      
out_res<-NULL        
for ( i in 1:length(genes)){
  
  gene_focal<-gene_pos[gene_pos$name==genes[i],]
  overlap<-aa[(aa$chrom == gene_focal$chrom & aa$pos >= gene_focal$start & aa$pos <= gene_focal$end),]
  overlap_count<-data.frame(rbind(table(overlap$CLS)))
  tt<-data.frame(t(overlap_count))
  
  tt$varint_type<-rownames(tt)
  
  tt$gene<-genes[i]
  out_res<-rbind.fill(tt,out_res)  
}

colnames(out_res)<-c("count","varint_type","gene")
write.table(out_res,file = "out_res1_driver_mutations_clonal_subclonal_cll.txt", col.names = TRUE, row.names = FALSE, quote= FALSE, sep = "\t")


ggplot(out_res, aes(fill=varint_type, y=count, x=gene)) + 
  geom_bar(position='stack', stat='identity') +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 9.5))



