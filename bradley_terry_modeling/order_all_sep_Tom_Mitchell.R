args=commandArgs(TRUE)
run<-as.integer(args[1])
input.folder <- "/lustre/scratch112/sanger/casm/cgp/tjm/Prostate/TCGA/Jan_2016"
output.folder  <- "/lustre/scratch112/sanger/casm/cgp/tjm/Prostate/TCGA/Jan_2016/results"

library("BradleyTerry2")
library("qvcalc")
library("gtools")

maxchr <- read.table(paste(input.folder,"/chrlengths.txt", sep=""))
maxchr <- maxchr[[1]]
loci <- read.table(paste(input.folder,"/cumchrlenths.txt", sep=""))
loci <- loci[[1]]
allsegsa <- read.table(paste(input.folder,"/allsegments.txt", sep=""), header=T)
allsegsa <- allsegsa[!allsegsa$chr == "X",]
allsegsa$chr <- as.numeric(allsegsa$chr)
#allsegsa$chr <- as.numeric(levels(allsegsa$chr))[allsegsa$chr]
t.e.l <-  c(39700000, 43000000)

ETS_pos <- read.table("/lustre/scratch112/sanger/casm/cgp/tjm/Prostate/TCGA/Jan_2016/TCGA_ETS+.txt")[[1]]
ETS_pos <- unique(substr(ETS_pos,9,12))
short_names <- substr(sub("X*","",allsegsa[,1]),1,4)

# set correct names for the run
status <- c("all","TMPRSS_positive", "TMPRSS_negative")
if (run==1) {names <- unique(allsegsa$Tumour_Name)}
if (run==2) {names <- unique(allsegsa$Tumour_Name[short_names%in%ETS_pos])}
if (run==3) {names<- unique(allsegsa$Tumour_Name[!short_names%in%ETS_pos])}



gain.segs <- allsegsa[grep("*Gain",allsegsa$CNA),]
remove.subclonal.gains <- gain.segs[duplicated(gain.segs[1:4]),]
allsegsa<-allsegsa[!row.names(allsegsa)%in%row.names(remove.subclonal.gains),]
loss.segs <- allsegsa[grep("*Loss",allsegsa$CNA),]
remove.subclonal.losses <- loss.segs[duplicated(loss.segs[1:4]),]
allsegsa<-allsegsa[!row.names(allsegsa)%in%row.names(remove.subclonal.losses),]

sigLOH <- read.table(paste(input.folder,"/CN_landscape/",status[run],"_samples_significant_pvals_FDR_LOH.txt",sep=""))
sigGain <- read.table(paste(input.folder,"/CN_landscape/",status[run],"_samples_significant_pvals_FDR_Gain.txt", sep=""))
sigHD <- read.table(paste(input.folder,"/CN_landscape/",status[run],"_samples_significant_pvals_FDR_HD.txt", sep=""))
sigLOH <- sigLOH[,c(1:4)]
sigGain <- sigGain[,1:4]
sigHD <- sigHD[,1:4]


# reduce segment length and number to >minrefseglength bp
CNAs <- c("sigLOH", "sigGain","sigHD")
out <- c("sigLOHc","sigGainc","sigHDc")
minrefseglength <- 1000000
for (k in 1:length(CNAs)) {
  cut <- NULL
  data <- get(CNAs[k])
  c <- unique(data$chr)
  for (j in 1:length(c)) {
    data.c <- data[data$chr == c[j],]
    if (dim(data.c)[1]>1) {
      for (i in 1:(dim(data.c)[1]-1)) {
        cuti <- data.c[i,]
        if ((cuti$ep-cuti$sp)>minrefseglength|(data.c$sp[i+1]-data.c$ep[i])>2000000) {
          cut <- rbind(cut, cuti)
        }
        else {
          data.c$sp[i+1] <- data.c$sp[i]
        }
      }
      cuti <- data.c[i+1,]
      cut <- rbind(cut, cuti)
    }
  }
  assign(out[k],cut)  
}

namL <- paste("L",sigLOHc$chr, ":", sigLOHc$sp/1000000,"-", sigLOHc$ep/1000000, sep="")
namG <- paste("G",sigGainc$chr, ":", sigGainc$sp/1000000,"-", sigGainc$ep/1000000, sep="")
namHD <- paste("HD",sigHDc$chr, ":", sigHDc$sp/1000000,"-", sigHDc$ep/1000000, sep="")
ln <- length(namL)+length(namG)+length(namHD)

clusterfiles = list.files(paste(input.folder,"/DP/clusters", sep=""), pattern="*Assignments.txt")
clustersummaries = list.files(path = paste(input.folder,"/DP/clusters", sep=""), pattern="*Positions.txt")
ordersamples <- sub("_withSubclonality_muts_withDPAssignments.txt","",clusterfiles)
scores <- matrix(rep(0,(ln+2)^2), nrow = ln+2, ncol = ln+2)
corr <- matrix(rep(0,(ln+2)^2), nrow = ln+2, ncol = ln+2)
dimnames(scores) <- list(c(namL, namG, namHD,"TMPRSS", "WGD"), c(namL, namG, namHD,"TMPRSS", "WGD"))
dimnames(corr) <- list(c(namL, namG, namHD,"TMPRSS", "WGD"), c(namL, namG, namHD,"TMPRSS", "WGD"))

for (i in 1:length(names)) {
  file <- paste(input.folder,"/DP/clusters/",names[i],"_withSubclonality_muts_withDPAssignments.txt", sep="")
  if (file.exists(file)) {
    subclonal <- read.table(file = file, header=T)
    add <- subclonal[subclonal$nMaj1_A==0&subclonal$nMin1_A==0,]
    if (dim(add)[1]>0) {
      add$CNA<-"sHD"
      subclonal <- rbind(subclonal, add)
    }
    summary <- read.table(file = paste(input.folder,"/DP/clusters/",names[i],"_withSubclonality_bestNodePositions.txt", sep=""), header=T)
    clonal <- allsegsa[allsegsa$Tumour_Name%in%names[i]&allsegsa$frac1_A==1,]
    fraction <- summary$subclonal.fraction[subclonal$DirichletProcessCluster]
    subclonal$frac1_A <- fraction
    subclonal$DirichletProcessCluster <- NULL
    allsegsi <- rbind(clonal, subclonal)
  }
  if (!file.exists(file)) {
    allsegsi <- allsegsa[allsegsa$Tumour_Name%in%names[i],]
  }
  
  gain.segs <- allsegsi[grep("*Gain",allsegsi$CNA),]
  remove.subclonal.gains <- gain.segs[duplicated(gain.segs[1:4]),]
  allsegsi<-allsegsi[!row.names(allsegsi)%in%row.names(remove.subclonal.gains),]
  loss.segs <- allsegsi[grep("*Loss",allsegsi$CNA),]
  remove.subclonal.losses <- loss.segs[duplicated(loss.segs[1:4]),]
  allsegsi<-allsegsi[!row.names(allsegsi)%in%row.names(remove.subclonal.losses),]
  
  incL <- NULL
  for (w in 1:length(sigLOHc$sp)) {
    inci <- (allsegsi$startpos < sigLOHc$ep[w]) & (allsegsi$endpos > sigLOHc$sp[w]) & (allsegsi$chr == sigLOHc$chr[w]) &
      (allsegsi$CNA == "cLoss" | allsegsi$CNA == "sLoss")
    incL <- cbind(incL, inci)
  }
  dimnames(incL)[[2]][1:(length(namL))] <- namL
  allsegsinc <- cbind(allsegsi, incL)
  
  incG <- NULL
  for (w in 1:length(sigGainc$sp)) {
    inci <- (allsegsi$startpos < sigGainc$ep[w]) & (allsegsi$endpos > sigGainc$sp[w]) & (allsegsi$chr == sigGainc$chr[w])&
      (allsegsi$CNA == "cGain" | allsegsi$CNA == "sGain" | allsegsi$CNA == "cAmp" | allsegsi$CNA == "sAmp")
    incG <- cbind(incG, inci)
  }
  dimnames(incG)[[2]][1:(length(namG))] <- namG
  allsegsinc <- cbind(allsegsinc, incG)
  
  incHD <- NULL
  for (w in 1:length(sigHDc$sp)) {
    inci <- (allsegsi$startpos < sigHDc$ep[w]) & (allsegsi$endpos > sigHDc$sp[w]) & (allsegsi$chr == sigHDc$chr[w])&
      (allsegsi$CNA == "cHD" | allsegsi$CNA == "sHD")
    incHD <- cbind(incHD, inci)
  }
  dimnames(incHD)[[2]][1:(length(namHD))] <- namHD
  allsegsinc <- cbind(allsegsinc, incHD)
  
  incTMPRSS <- (allsegsi$startpos>t.e.l[1])&(allsegsi$endpos<t.e.l[2])&(allsegsi$chr==21)&
    (allsegsi$CNA == "cLOH"|allsegsi$CNA == "sLOH")
  allsegsinc <- cbind(allsegsinc, incTMPRSS)
  
  # filter out segments that do not cross reference segments and do not have CNVs of interest
  include <- apply(allsegsinc[,c(11:(11+ln))]==TRUE,1,any)
  allsegsinc <- allsegsinc[include,]
  allsegsinc <- allsegsinc[allsegsinc$CNA!="NoCNV"&allsegsinc$CNA!="NotDet",]
  
  corri <- matrix(0, nrow = ln+2, ncol = ln+2)
  scoresi <- matrix(0, nrow = ln+2, ncol = ln+2)
  ordinta = NULL
  ordintna = NULL
  za <- NULL
  maxin <- 0
  ordi <- allsegsinc
  # This allocates a point to scores if a given event occurs prior to another event using the pigeonhole principle for clonal versus subclonal events
  if (dim(ordi)[1] > 1){
    # this is the subclonal ordering script
    ordi <- ordi[order(ordi$frac1_A, decreasing = T),]
    ch_all <- ordi[,c(11:(11+ln))]
    for (j in 2:dim(ordi)[1]) {
      maxin <- c(maxin, max(c((ordi$frac1_A[j-1]-ordi$frac1_A[j]),maxin)))
    }
    if (ordi$tumour_ploidy[1]>2) {
      corri[ln+2,ln+2] <- corri[ln+2,ln+2]+1
    }
    for (n in 1:dim(ch_all)[1]) {
      for (p in 1:dim(ch_all)[1]) {
        xfrac <- ordi$frac1_A[p]
        yfrac <- ordi$frac1_A[n]
        xgap <- maxin[p]
        nc <- as.logical(c(ch_all[n,1:(ln+1)],FALSE))
        pc <- as.logical(c(ch_all[p,1:(ln+1)],FALSE))
        corri[nc,pc] <- corri[nc,pc]+1
        if ((xfrac>yfrac) & (xgap<yfrac)) {
          scoresi[nc,pc] <- scoresi[nc,pc] +1
        }
      }
    }
    # new:  if there is clonal HD, there is an implicit win against HD for loss. We do not know at present which of the two LOH hits occurred first
    # and which resulting in the HD segment.
    l.HD <- which(ordi$CNA=="cHD")
    if (length(l.HD)>0) {
      for (l in 1:length(l.HD)) {
        pc<-as.logical(c(ch_all[l.HD[l]+1,],FALSE))
        nc<-as.logical(c(ch_all[l.HD[l],],FALSE))
        scoresi[nc,pc] <- scoresi[nc,pc] +1
      }
    }
    # new: time WGD
    WGD.include <- ordi[ordi$tumour_ploidy==4,]
    if (dim(WGD.include)[1]>0) {
      ch_all <- ordi[,c(11:(11+ln))]
      for (j in 1:dim(WGD.include)[1]) {
        pc <- as.logical(c(ch_all[j,1:(ln+1)],FALSE))
        nc <- as.logical(c(ch_all[j,1:(ln+1)],FALSE))
        corri[ln+2,pc] <- corri[ln+2,pc] +1            
        corri[nc,ln+2] <- corri[nc,ln+2] +1
        # loss/ gain first then WGD
        if (WGD.include$CNA[j]=="cLoss"  & WGD.include$nMin1_A[j]==0) {
          scoresi[ln+2,pc] <- scoresi[ln+2,pc] +1            
        }
        if (WGD.include$CNA[j]=="cLoss"  & WGD.include$nMin1_A[j]==1) {
          scoresi[nc,ln+2] <- scoresi[nc,ln+2] +1
        }
        if (WGD.include$CNA[j]=="cHD") {
          scoresi[ln+2,pc] <- scoresi[ln+2,pc] +1
        }
        if (WGD.include$CNA[j]=="cGain" & WGD.include$nMaj1_A[j]>3) {
          scoresi[ln+2,pc] <- scoresi[ln+2,pc] +1
        }
        if (WGD.include$CNA[j]=="cGain" & WGD.include$nMaj1_A[j]==3) {
          scoresi[nc,ln+2] <- scoresi[nc,ln+2] +1
        }
      }        
    }
  }
  if (sum(grepl(substr(sub("X*","",names[i]),1,4),ETS_pos))>0) {
    corri[ln+1,] <- (colSums(corri)>0) +0
    corri[,ln+1] <- (colSums(corri)>0) +0
    corri[ln+1,ln+1] <- 1
  }
  corr[corri>0]<-corr[corri>0]+1
  scores[scoresi>0]<-scores[scoresi>0]+1  
  print(names[i])
}

sigc <- rbind(cbind("CNA"="LOH", sigLOHc), cbind("CNA"="Gain", sigGainc), cbind("CNA"="HD", sigHDc),cbind("CNA"="TMPRSS-ERG", "chr"=21, "sp"=t.e.l[1], "ep"=t.e.l[2], "val"=corr[ln+1,ln+1]),
              cbind("CNA"="WGD", "chr"=NA, "sp"=NA, "ep"=NA, "val"=length(unique(allsegsa[allsegsa$tumour_ploidy>2,1]))))
colnames(scores) <- c(namL,namG, namHD,"TMPRSS2-ERG", "WGD")
write.table(scores, file = paste(output.folder , "/scores_",status[run],"_sep.txt", sep=""),sep="\t")
write.table(corr, file = paste(output.folder , "/corr_",status[run],"_sep.txt", sep=""),sep="\t")



sing<-which(rowSums(scores)>00|colSums(scores)>00)
scores<- scores[sing,sing]
sigc <- sigc[sing,]
for (i in 1:dim(corr)[1]) {
  corr[i,]<-corr[i,]/corr[i,i]
}

players <- rownames(scores)
options(expressions = 500000)
comb <- combinations(nrow(scores), 2)
won <- scores[comb]
lost <- t(scores)[comb]
res <- !(won == 0 & lost == 0)
player1 <- factor(players[comb[, 1]], levels = players)[res]
player2 <- factor(players[comb[, 2]], levels = players)[res]
scores.sf <- data.frame(player1, player2, win1 = won[res], win2 = lost[res])
names(scores.sf)[1:2] <- c("loci_1", "loci_2")

orderModel <- BTm(cbind(win1, win2), loci_1, loci_2, ~ loci_,
                  id = "loci_", data = scores.sf, br = T, na.action = na.omit)
print(summary(orderModel))
save(orderModel, file = paste(output.folder , "/orderModel_",status[run],"_sep.Rdata", sep=""))

orderModel.qv <- BTabilities(orderModel)

orderModel.qv <- cbind(orderModel.qv, sigc[1:dim(sigc)[1],])
write.table(orderModel.qv, file = paste(output.folder , "/orderModel_",status[run],"_sep.txt", sep=""), sep="\t")
