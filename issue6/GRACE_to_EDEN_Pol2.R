## GRACE to EDEN plots ##
## Berkley Gryder and  ##
## Hyunmin Kim,  2022  ##
## Bhava, 2023         ##
#Boxplot of Log2FC in gene expression before and after dCBP1, for genes near p300 enhancers with and withour P3F

#working directory
setwd("/Volumes/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2_real/GRACE_EDEN_Pol2/")

#load libraries
library(plyr);library(ggplot2);library(dplyr)

#Set Variables
Sample_Name<- "RH4_Pol2_GRACE"

#EDEN.GRACE
##Data from EDEN_multi_genes.txt output
EDEN.near <- read.table("RH4_DMSO_dCBP1_merged_Pol28k_3m_2.wP3F_GRACE_TPM5_nearest-genes.txt", sep = "\t", header = T)
EDEN.over <- read.table("RH4_DMSO_dCBP1_merged_Pol28k_3m_2.wP3F_GRACE_TPM5_overlap-genes.txt", sep = "\t", header = T)
colnames(EDEN.near)=c("chr","start","end","shortAQuA","longAQuA","region.name","loops","group","rank.CPM","wP3F","Gene","TPM","distance")
colnames(EDEN.over)=c("chr","start","end","shortAQuA","longAQuA","region.name","loops","group","rank.CPM","wP3F","Gene","distance","TPM")
EDEN.over=EDEN.over[,c("chr","start","end","shortAQuA","longAQuA","region.name","loops","group","rank.CPM","wP3F","Gene","TPM","distance")]

EDEN.fill = subset(EDEN.near, !(EDEN.near$region.name %in% EDEN.over$region.name))
EDEN.full = rbind(EDEN.over,EDEN.fill)

##remove rows with no associated gene name
EDEN.GRACE <- subset(EDEN.full, Gene!="."); EDEN.GRACE$TPM = as.numeric(EDEN.GRACE$TPM)

##Join the gene expression matrix
EXP = read.table("RH4_RH5_p300.matrix.txt", header = T)

EDEN.GRACE.EXP = join(EDEN.GRACE,EXP,by = 'Gene') 


##Read the new bed file that has both p300 and P3F binding information
##EDEN.GRACE p300
EDEN.GRACE.p300<-read.table("RH4_Pol2_GRACE.EDEN.p300_peaks.bed")

colnames(EDEN.GRACE.p300)=c("chr","start","end","wP3F","wp300")

##making a new coordinates column 
EDEN.GRACE.p300$coord = paste0(EDEN.GRACE.p300$chr,":",EDEN.GRACE.p300$start,"-",EDEN.GRACE.p300$end)

EDEN.GRACE$coord<-paste0(EDEN.GRACE$chr,":",EDEN.GRACE$start,"-",EDEN.GRACE$end)

##using coordinates to join the new bed file and the EDEN.GRACE file used above that has the TPM, log2fc and AQuA values
EDEN.GRACE.complete<-join(EDEN.GRACE.p300,EDEN.GRACE, by = 'coord') 

##6091 rows were duplicated
EDEN.GRACE.complete <- EDEN.GRACE.complete[!duplicated(EDEN.GRACE.complete), ]

##keeping only the columns needed downstream
EDEN.GRACE.complete<-EDEN.GRACE.complete[,c(1,2,3,4,5,6,10,11,12,13,14,15,17,18,19)]

##Join RH4/RH5_DMSO/dCBP1 gene expression matrix to the EDEN.GRACE matrix
EDEN.GRACE.EXP.groups = join(EDEN.GRACE.complete,EXP,by = 'Gene') 
# ntile = 10
# EDEN.GRACE.EXP.groups$Pol2CPM_ntile <-  as.character( sprintf("%02d", ntile(EDEN.GRACE.EXP.groups$sumAQuA, ntile))  )

##adding a column that groups data into p300 w/o P3F and p300 w P3F based on the values in the p300 and P3F columns
EDEN.GRACE.EXP.groups$group_p300.P3F <- ifelse(EDEN.GRACE.EXP.groups$wp300 > 0 & EDEN.GRACE.EXP.groups$wP3F == 0, "p300 without P3F", "p300 with P3F")
colnames(EDEN.GRACE.EXP.groups)

# colnames(EDEN.GRACE.EXP.groups)<-c("chr", "start", "end", "wP3F", "wp300", "coord", "shortAQuA", "longAQuA", "sumAQuA", 
#                                    "region.name", "loops", "group", "rank.CPM", "Gene", "TPM", "distance", 
#                                    "RH4_A485_072221_CWRU", "RH4_dCBP1_072221_CWRU", "RH4_DMSO_072221_CWRU", "RH4_NT_072221_CWRU", 
#                                    "RH5_A485_072221_CWRU", "RH5_dCBP1_072221_CWRU", "RH5_DMSO_072221_CWRU",  "RH5_NT_072221_CWRU", 
#                                    "L2FC_RH4_A485", "L2FC_RH4_dCBP1", "L2FC_RH4_NT", "L2FC_RH5_A485", "L2FC_RH5_dCBP1", "L2FC_RH5_NT", "L2FC_A485", "L2FC_dCBP1", "L2FC_p300", 
#                                    "PVAL_A485", "PVAL_dCBP1", "Pol2CPM_ntile", "group_p300.P3F")

#plot it
require(graphics);library(ggplot2);library(RColorBrewer);library(plyr);library(dplyr)
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
mycols = colorRampPalette(my3cols)(ntile)

class(EDEN.GRACE.EXP.groups$L2FC_RH4_dCPB1)
class(EDEN.GRACE.EXP.groups$group_p300.P3F)
#RH4 - dCBP1
#ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=log2(TPM+1),color=group))+geom_jitter()+scale_color_manual(values = my3cols)+theme_bw()+geom_boxplot(outlier.shape = NA)
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH4_dCPB1,group=group_p300.P3F,fill=group))+geom_jitter(alpha=0.05,aes(color=group))+
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

#Just merged clusters

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

#############################################################################
#Grouping based on genes rather than 3D information

housekeeping.genes<-read.delim("HouseKeeping_genes.txt", sep = "\t", header = F)
housekeeping.genes<-housekeeping.genes$V1

crc.genes<-read.delim("GRYDER_RH4_CRC_TFs_TOP50prct.genelist.txt", sep = "\t", header = F)
crc.genes<-crc.genes$V1

EDEN.GRACE.EXP.groups <- EDEN.GRACE.EXP.groups %>%
  mutate(group_genes = case_when(
    Gene %in% housekeeping.genes ~ "Housekeeping genes",
    Gene %in% crc.genes ~ "CRC genes",
    TRUE ~ "Other p300 binding genes"
  ))

my2cols <- c("#007A86", "#25B1BB")
mycols = colorRampPalette(my2cols)(ntile)

##Based on the gene groups
##RH4 - dCBP1
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH4_dCPB1,group=group_p300.P3F,fill=group_p300.P3F))+geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my2cols)+scale_color_manual(values = my2cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group_genes)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

#RH4 - A485
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH4_A485,group=group_p300.P3F, fill=group_p300.P3F))+geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my2cols)+scale_color_manual(values = my2cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group_genes)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

##RH5 - dCBP1
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH5_dCBP1,group=group_p300.P3F,fill=group_p300.P3F))+geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my2cols)+scale_color_manual(values = my2cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group_genes)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

#RH5 - A485
ggplot(EDEN.GRACE.EXP.groups,aes(x=group_p300.P3F, y=L2FC_RH5_A485,group=group_p300.P3F,fill=group_p300.P3F))+geom_jitter(alpha=0.05,aes(color=group_p300.P3F))+
  scale_fill_manual(values = my2cols)+scale_color_manual(values = my2cols)+theme_bw()+geom_boxplot(outlier.shape = NA)+facet_wrap(~group_genes)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

#################################################################################
