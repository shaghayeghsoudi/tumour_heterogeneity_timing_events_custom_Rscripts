#!/usr/bin/env Rscript

#rm(list = ls())
## this script runs TrackSig package to track changes in muttaional influences during the process of tumour evolution in a parallel mode in R.


##### ATTRIBUTION #####
# Original Author:  Shaghayegh Soudi
# Contributors:    NA 


library(foreach)
library(doParallel)
#detectCores()
library(stringr)
library("tidyverse")


## path to Sigprofiler results to take the first four/five active signatures for each patient
## Tracksig has a function that detects active signature but the developer suggests to do it with a more proper package. We ran SigProfiler fisrt to find active signatures per patient and used the output for the TrackSig nalysis

cosmic<-read.table(file = "/PATH/TO/SIGPROFILER/OUTPUTS/COSMIC_SBS96_Activities_refit.txt", header = TRUE)
cosmic_with_rname<-cosmic %>% remove_rownames %>% column_to_rownames(var="Samples") 


shared_samples<-read.table(file = "shared_samples.table", header = TRUE)
tumour<-data.frame(str_split_fixed(shared_samples$shared_samples,"--",2))
colnames(tumour)<-c("tumour","normal")

shared_samples<-cbind(shared_samples,tumour)
shared_samples$tumour_type<-rep("CLL")
#shared_samples<-as.character(shared_samples$shared_samples)

### directories of merged vcfs and CNVs 
vcf_dir<-"/PATH/TO/vcfs"
cna_dir<-"/PATH/TO/cna"

contamination<-read.table(file = "sample_normal_contamination.txt", header = TRUE)


shared_samples_test<-shared_samples[1:3,]

parallel::detectCores()
n.cores <- parallel::detectCores() - 50

### **create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

## **register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()


#for (j in 1:nrow(shared_samples)){


foreach(j = 1:nrow(shared_samples_test))%dopar%{   
  library(TrackSig)
  library(ggplot2)
  library(vcfR)
  library(data.table)
  #library(stringr)
  
  focal_shared_samples<-shared_samples[j,]
  
  
  focal_vcf<-paste(paste(vcf_dir,focal_shared_samples$shared_samples, sep = "/"),"vcf",sep = ".")
  vcf <- read.vcfR(focal_vcf, verbose = FALSE )
  
  if(nrow(vcf)> 500){
    
    
    
    focal_cna<-paste(paste(cna_dir,focal_shared_samples$shared_samples, sep = "/"),"_cna.txt",sep = "")
    contamination_focal<-contamination[contamination$SAMPLE==as.character(focal_shared_samples$shared_samples),2]
    
    
    #detectedSigs <- detectActiveSignatures(vcfFile = focal_vcf, cnaFile = focal_cna, purity = contamination_focal, threshold = 0.05)
    #detectedSigs <- detectActiveSignatures(vcfFile = focal_vcf, cnaFile = focal_cna, purity = contamination_focal, threshold = 0.05)
    
    
    focal_sig<-cosmic_with_rname[rownames(cosmic_with_rname)==as.character(focal_shared_samples$tumour),]
    focal_sig_active<-focal_sig[,which(focal_sig[1,]>0)]
    
    if (ncol(focal_sig_active) <= 4){
      focal_sig_active_selected<-focal_sig[,which(focal_sig[1,]>0)] 
    }
    
    if (ncol(focal_sig_active)> 4){
      focal_sig_active_selected<-sort(focal_sig_active, decreasing = TRUE)[,1:4]
      
    }
    
    
    detectedSigs <- colnames(focal_sig_active_selected)
    
    binSize = c("100","200","300")
    
    for (bb in 1:length(binSize)){
      
      scoreMethod = "SigFreq"
      
      
      
      traj <- TrackSig(sampleID = as.character(focal_shared_samples$shared_samples), activeInSample = detectedSigs,
                       vcfFile = focal_vcf, cnaFile = focal_cna, purity = contamination_focal,
                       scoreMethod = scoreMethod, binSize = as.numeric(binSize[bb]))
      
      
      write.table(traj$mixtures, file = paste("traj_mixture_",shared_samples[j,1],"_bin_size",binSize[bb],".table", sep = ""), col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
      write.table(traj$binData, file = paste("traj_binData_",shared_samples[j,1],"_bin_size",binSize[bb],".table", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      
      change_point<-data.frame(traj$changepoints)
      write.table(change_point, file = paste("change_points",shared_samples[j,1],"_bin_size",binSize[bb],".table", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      
      
      ## plot 
      VAF_coding_plot= paste("traj_linear_x-axis",shared_samples[j,1],"_bin_size",binSize[bb],".pdf", sep = "")
      pdf(VAF_coding_plot, width = 12, height =4)
      plotTrajectory<-plotTrajectory(traj, linearX = T) + labs(title = paste(focal_shared_samples[,4],focal_shared_samples[,1]), sep = "")
      grid::grid.draw(plotTrajectory)
      dev.off() 
      
      
      VAF_coding_plot_nonl= paste("traj_non-linear_x-axis",shared_samples[j,1],"_bin_size",binSize[bb],"_nonLinPlot.pdf", sep = "")
      pdf(VAF_coding_plot_nonl, width = 12, height =4)
      nonLinPlot<-plotTrajectory(traj, linearX = F, anmac = T) + labs(paste(focal_shared_samples[,4],focal_shared_samples[,1]), sep = "")
      grid::grid.draw(nonLinPlot)
      dev.off()
      
      
      VAF_coding_plot_hist= paste("traj_",shared_samples[j,1],"_bin_size",binSize[bb],"_hist.pdf", sep = "")
      pdf(VAF_coding_plot_hist, width = 12, height =6)
      plotGrobs <- addPhiHist(traj, nonLinPlot)
      grid::grid.draw(plotGrobs)
      dev.off()
      
      
     
      
    } ## bin
    
  } # if
  
} # j

parallel::stopCluster(cl = my.cluster)


