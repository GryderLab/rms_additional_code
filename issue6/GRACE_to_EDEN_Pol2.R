## GRACE to EDEN plots ##
## Berkley Gryder and  ##
## Hyunmin Kim,  2022  ##
#########################

## Read in GRACE EDEN outputs
## Plotting data from EDEN ##
library(plyr);library(ggplot2);library(dplyr)

setwd("/Volumes/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2_real/GRACE_EDEN_Pol2/")
#Set Variables
Sample_Name<- "RH4_Pol2_GRACE"

## read in EDEN _multi-genes.txt output file
EDEN.near <- read.table("RH4_DMSO_dCBP1_merged_Pol28k_3m_2.wP3F_GRACE_TPM5_nearest-genes.txt", sep = "\t", header = T)
EDEN.over <- read.table("RH4_DMSO_dCBP1_merged_Pol28k_3m_2.wP3F_GRACE_TPM5_overlap-genes.txt", sep = "\t", header = T)
colnames(EDEN.near)=c("chr","start","end","shortAQuA","longAQuA","region.name","loops","group","rank.CPM","wP3F","Gene","TPM","distance")
colnames(EDEN.over)=c("chr","start","end","shortAQuA","longAQuA","region.name","loops","group","rank.CPM","wP3F","Gene","distance","TPM")
EDEN.over=EDEN.over[,c("chr","start","end","shortAQuA","longAQuA","region.name","loops","group","rank.CPM","wP3F","Gene","TPM","distance")]

EDEN.fill = subset(EDEN.near, !(EDEN.near$region.name %in% EDEN.over$region.name))
EDEN.full = rbind(EDEN.over,EDEN.fill)

## remove rows with no associated gene name
EDEN.GRACE <- subset(EDEN.full, Gene!="."); EDEN.GRACE$TPM = as.numeric(EDEN.GRACE$TPM)
EDEN.GRACE$sumAQuA = EDEN.GRACE$shortAQuA + EDEN.GRACE$longAQuA

## Create bed; check peak location in IGV - columns = chromosome, start, end, group
write.table(EDEN.GRACE[,c(1,2,3,10)],file = paste(Sample_Name,".EDEN.bed",sep = ""), col.names = F, row.names = F,quote = F)

## join a gene expression matrix
EXP = read.table("RH4_RH5_p300.matrix.txt", header = T)
#EXP = EXP[,c(grep("Gene|LNCaP",names(EXP))), drop=F];names(EXP)[names(EXP) == 'GeneID'] <- 'Gene'

## Ntile classify GRACE data
EDEN.GRACE.EXP = join(EDEN.GRACE,EXP,by = 'Gene') 
ntile = 10
EDEN.GRACE.EXP$Pol2CPM_ntile <-  as.character( sprintf("%02d", ntile(EDEN.GRACE.EXP$sumAQuA, ntile))  )

#plot it
require(graphics);library(ggplot2);library(RColorBrewer);library(plyr);library(dplyr)
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
mycols = colorRampPalette(my3cols)(ntile)

ggplot(EDEN.GRACE.EXP,aes(x=group, y=log2(TPM+1),color=group))+geom_jitter()+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)#+facet_wrap(~group)
ggplot(EDEN.GRACE.EXP,aes(x=wP3F, y=L2FC_RH4_dCPB1,group=wP3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)

ggplot(EDEN.GRACE.EXP,aes(x=Pol2CPM_ntile, y=log2(TPM+1),color=Pol2CPM_ntile))+geom_point(alpha=0.1)+scale_color_manual(values = mycols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)
ggplot(EDEN.GRACE.EXP,aes(x=wP3F, y=log2(longAQuA+1),group=wP3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)

write.table(EDEN.GRACE.EXP,file = paste(Sample_Name,".EDEN.GRACE.EXP.txt",sep = ""), col.names = T, row.names = F,quote = F,sep = "\t")
######################################################################################
EDEN.GRACE.P3F = subset(EDEN.GRACE, EDEN.GRACE$wP3F > 0)
write.table(unique(EDEN.GRACE.P3F$Gene),file = paste(Sample_Name,".wP3F.genelist.txt",sep = ""), col.names = F, row.names = F,quote = F)

######################################################################################
## Check Gene overlapping (by Hyunmin)
i=rep(F,nrow(EDEN.GRACE.EXP)); i[grep("overlaps",EDEN.GRACE.EXP$group)]=T
g1=unique(EDEN.GRACE.EXP[i ,]$Gene );g2=unique(EDEN.GRACE.EXP[!i,]$Gene )
library(ggvenn)
ggvenn(list(overlap=g1,not=g2))
g12=intersect(g1,g2)

# add "overlap" genes from "not" category if the gene associates with any AR overlap peak (Berkley)
EDEN.GRACE.EXP$GeneOverlap = EDEN.GRACE.EXP$Gene %in% g1
ggplot(EDEN.GRACE.EXP,aes(x=as.character(H3K27ac_DeltaBG15n_5tile),y=TPM_L2FC_LNCaPXIP_BG15n,color=GeneOverlap))+geom_boxplot()
######################################################################################

## collapse EDEN.GRACE by Gene by Hyunmin
x=EDEN.GRACE.EXP%>%group_by(Gene)%>%
  dplyr::summarize( 
    Chromosome=Chromosome,Start=min(Start),End=max(End),
    group=paste(unique(sort(group)),collapse=","), 
    shortAQuA= sum(shortAQuA),
    longAQuA = sum(longAQuA),
    sumAQuA = sum(sumAQuA),
    TPM = as.numeric(TPM),
    rank.CPM = max(rank.CPM),
    TPM_L2FC_LNCaPXIP_BG15n = TPM_L2FC_LNCaPXIP_BG15n
  ) 

#x[grep(",",x$group),]$group = "both"
EDEN.GRACE.sum=unique(x)

ntile = 10
EDEN.GRACE.sum$GRACE_ntile <-  as.character( sprintf("%02d", ntile(EDEN.GRACE.sum$rank.CPM, ntile))  )

#plot it
require(graphics);library(ggplot2);library(RColorBrewer);library(plyr);library(dplyr)
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
mycolors = colorRampPalette(my3cols)(ntile)

ggplot(EDEN.GRACE.sum,aes(x=rank.CPM, y=log2(TPM+1),color=GRACE_ntile))+geom_point(alpha=0.5)+geom_boxplot(outlier.shape = NA)+scale_color_manual(values = mycolors)+theme_bw()
ggplot(EDEN.GRACE.sum,aes(x=rank.CPM, y=log2(TPM+1),color=GRACE_ntile))+geom_boxplot()+scale_color_manual(values = mycolors)+theme_bw()
ggplot(EDEN.GRACE.sum,aes(x=GRACE_ntile, y=TPM_L2FC_LNCaPXIP_BG15n,fill=GRACE_ntile))+geom_violin(width=0.8,alpha=0.5)+scale_fill_manual(values = mycolors)+theme_bw()+geom_boxplot(width=0.7,outlier.shape = NA)

##############################################################################
#including p300 data - Bhava
##Read the new bed file that has both p300 and P3F binding information
EDEN.GRACE.p300<-read.table("RH4_Pol2_GRACE.EDEN.p300_peaks.bed")

colnames(EDEN.GRACE.p300)=c("chr","start","end","wP3F","wp300")

##making a new coordinates column 
EDEN.GRACE.p300$coord = paste0(EDEN.GRACE.p300$chr,":",EDEN.GRACE.p300$start,"-",EDEN.GRACE.p300$end)

##making a new coordinates column 
EDEN.GRACE$coord<-paste0(EDEN.GRACE$chr,":",EDEN.GRACE$start,"-",EDEN.GRACE$end)

##using coordinates to join the new bed file and the EDEN.GRACE file used above that has the TPM, log2fc and AQuA values
EDEN.GRACE.complete<-join(EDEN.GRACE.p300,EDEN.GRACE, by = 'coord') 

##~6000 rows were duplicated
EDEN.GRACE.complete <- EDEN.GRACE.complete[!duplicated(EDEN.GRACE.complete), ]

##keeping only the columns needed downstream
EDEN.GRACE.complete<-EDEN.GRACE.complete[, c(1,2,3,4,5,6,10,11,20,12,13,14,15,17,18,19)]

##Ntile
EDEN.GRACE.EXP.groups = join(EDEN.GRACE.complete,EXP,by = 'Gene') 
ntile = 10
EDEN.GRACE.EXP.groups$Pol2CPM_ntile <-  as.character( sprintf("%02d", ntile(EDEN.GRACE.EXP.groups$sumAQuA, ntile))  )

##adding a column that groups data into p300 w/o P3F and p300 w P3F based on the values in the p300 and P3F columns
EDEN.GRACE.EXP.groups$group_p300.P3F <- ifelse(EDEN.GRACE.EXP.groups$wp300 > 0 & EDEN.GRACE.EXP.groups$wP3F == 0, "p300 without P3F", "p300 with P3F")

colnames(EDEN.GRACE.EXP.groups)<-c("chr", "start", "end", "wP3F", "wp300", "coord", "shortAQuA", "longAQuA", "sumAQuA", 
                                   "region.name", "loops", "group", "rank.CPM", "Gene", "TPM", "distance", 
                                   "RH4_A485_072221_CWRU", "RH4_dCBP1_072221_CWRU", "RH4_DMSO_072221_CWRU", "RH4_NT_072221_CWRU", 
                                   "RH5_A485_072221_CWRU", "RH5_dCBP1_072221_CWRU", "RH5_DMSO_072221_CWRU",  "RH5_NT_072221_CWRU", 
                                   "L2FC_RH4_A485", "L2FC_RH4_dCBP1", "L2FC_RH4_NT", "L2FC_RH5_A485", "L2FC_RH5_dCBP1", "L2FC_RH5_NT", "L2FC_A485", "L2FC_dCBP1", "L2FC_p300", 
                                   "PVAL_A485", "PVAL_dCBP1", "Pol2CPM_ntile", "group_p300.P3F")

#plot it
require(graphics);library(ggplot2);library(RColorBrewer);library(plyr);library(dplyr)
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
mycols = colorRampPalette(my3cols)(ntile)


#RH4 - dCBP1
#ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=log2(TPM+1),color=group))+geom_jitter()+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH4_dCBP1,group=group_p300.P3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#RH4 - A485
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH4_A485,group=group_p300.P3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##RH5 - dCBP1
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH5_dCBP1,group=group_p300.P3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#RH5 - A485
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH5_A485,group=group_p300.P3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#Did not use RH5_pol2_P3F/p300 input file

#without the groups
#RH4 - dCBP1
ggplot(EDEN.GRACE.EXP.groups, aes(x = group_p300.P3F, y = L2FC_RH4_dCBP1, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)
#RH4 - A485
ggplot(EDEN.GRACE.EXP.groups, aes(x = group_p300.P3F, y = L2FC_RH4_A485, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)
##RH5 - dCBP1
ggplot(EDEN.GRACE.EXP.groups, aes(x = group_p300.P3F, y = L2FC_RH5_dCBP1, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)
#RH5 - A485
ggplot(EDEN.GRACE.EXP.groups, aes(x = group_p300.P3F, y = L2FC_RH5_A485, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)


#merged clusters

EDEN.GRACE.EXP.merged.culsters<-subset(EDEN.GRACE.EXP.groups, group == "merged clusters =1511")

#RH4 - dCBP1
ggplot(EDEN.GRACE.EXP.merged.culsters, aes(x = group_p300.P3F, y = L2FC_RH4_dCBP1, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+
  ggtitle("RH4_dCBP1: L2FC values - merged clusters")
#RH4 - A485
ggplot(EDEN.GRACE.EXP.merged.culsters, aes(x = group_p300.P3F, y = L2FC_RH4_A485, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+
  ggtitle("RH4_A485: L2FC values - merged clusters")
##RH5 - dCBP1
ggplot(EDEN.GRACE.EXP.merged.culsters, aes(x = group_p300.P3F, y = L2FC_RH5_dCBP1, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+
  ggtitle("RH5_dCBP1: L2FC values - merged clusters")
#RH5 - A485
ggplot(EDEN.GRACE.EXP.merged.culsters, aes(x = group_p300.P3F, y = L2FC_RH5_A485, fill = group_p300.P3F)) +
  geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my3cols)+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+
  ggtitle("RH5_A485: L2FC values - merged clusters")



################################################################################
# ##new bed file - 3 groups
# EDEN.GRACE_grp<-read.table("RH4_3groups.bed")
# 
# colnames(EDEN.GRACE_grp)=c("chr","start","end","wP3F","wp300")
# 
# ##making a new coordinates column 
# EDEN.GRACE_grp$coord = paste0(EDEN.GRACE_grp$chr,":",EDEN.GRACE_grp$start,"-",EDEN.GRACE_grp$end)
# 
# ##making a new coordinates column 
# EDEN.GRACE$coord<-paste0(EDEN.GRACE$chr,":",EDEN.GRACE$start,"-",EDEN.GRACE$end)
# 
# ##using coordinates to join the new bed file and the EDEN.GRACE file used above that has the TPM, log2fc and AQuA values
# EDEN.GRACE.complete_grp<-join(EDEN.GRACE_grp,EDEN.GRACE, by = 'coord') 
##EDEN.GRACE doesn't have the coordinated for no_P3F/p300 groups


