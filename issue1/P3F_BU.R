getwd()
setwd("/Volumes/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2_real/P3F_motifs")

# setwd("/Volumes/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2/P3F_motifs")
##lists of genes that are dependent or independent on C793 (from the paper) 
###two gene lists

C793_dependent<-read.delim("C793S_dependent_genes.txt", header = F, sep = " ")
C793_dependent<-C793_dependent$V1
C793_low.dependent<-read.delim("C793S_partial_rescue_independent_genes.txt", header = F, sep = " ")
C793_low.dependent<-C793_low.dependent$V1

#get the gene info from bed file for the genes
##bed file
Sample_Name = "P3F"
Folder_Name = "P3F_motifs"
GENE_bed<-read.table("gene_tss5k.clean.bed")
colnames(GENE_bed)<-c("Chrom", "Start", "End", "Gene", "Score", "Strand")

#filtering for the genes in the genelist
##C793_dependent
C793_dependent.genetable <- GENE_bed[GENE_bed$Gene %in% C793_dependent, ]
missing_genes_1 <- C793_dependent[!(C793_dependent %in% C793_dependent.genetable$Gene)]

##way to identify repeats
# library(dplyr)
# 
# # Identify genes repeated twice
# repeated_genes <- C793_dependent.genetable %>%
#   group_by(Gene) %>%
#   filter(n() == 2) %>%
#   distinct(Gene)
# repeated_genes$Gene

#"RIMBP3B" gene repeated twice

write.table(C793_dependent.genetable, "TSS_C793-dependent.genelist.bed",col.names = F, row.names = F, quote = F, sep="\t")

##C793_low.dependent
C793_low.dependent.genetable<-GENE_bed[GENE_bed$Gene %in% C793_low.dependent, ]
missing_genes_2 <- C793_low.dependent[!(C793_low.dependent %in% C793_low.dependent.genetable$Gene)]

# # Identify genes repeated twice
# repeated_genes <- C793_low.dependent.genetable %>%
#   group_by(Gene) %>%
#   filter(n() == 2) %>%
#   distinct(Gene)
# repeated_genes$Gene

#ICOSLG gene repeated twice

write.table(C793_low.dependent.genetable, "TSS_C793-low.dependent.genelist.bed",col.names = F, row.names = F, quote = F, sep="\t")

###############################################################################
##filtering for genes within 500kb from the bedtools intersect output
##C793 low dependent
P3F_C793_low.dependent<-read.table("P3F_near_C793_low.dependent.bed")

##rename bedtools closest col names
# colnames(P3F_C793dependent) <- c("chrom_A", "start_A", "end_A", "name_A", "score_A", "strand_A", "chrom_B", "start_B", "end_B", "name_B", "score_B", "strand_B", "distance")

req_distance=500000
P3F_C793_low.dependent$distance=P3F_C793_low.dependent[,ncol(P3F_C793_low.dependent)]

P3F_C793_low.dep_data <- subset(P3F_C793_low.dependent, P3F_C793_low.dependent$distance >= -(req_distance) & P3F_C793_low.dependent$distance <= req_distance)

P3F_C793_low.dep_data <- subset(P3F_C793_low.dep_data, P3F_C793_low.dep_data$distance != -1)
write.table(P3F_C793_low.dep_data, "P3F_C793_low.dep_data.bed",col.names = F, row.names = F, quote = F, sep="\t")

P3F_C793_low.dependent_genelist<-unique(P3F_C793_low.dep_data$V14)

##C793 dependent
P3F_C793_dependent<-read.table("P3F_near_C793dependent.bed")

##rename bedtools closest col names
# colnames(P3F_C793dependent) <- c("chrom_A", "start_A", "end_A", "name_A", "score_A", "strand_A", "chrom_B", "start_B", "end_B", "name_B", "score_B", "strand_B", "distance")

P3F_C793_dependent$distance=P3F_C793_dependent[,ncol(P3F_C793_dependent)]

P3F_C793_dep_data <- subset(P3F_C793_dependent, P3F_C793_dependent$distance >= -(req_distance) & P3F_C793_dependent$distance <= req_distance)

P3F_C793_dep_data <- subset(P3F_C793_dep_data, P3F_C793_dep_data$distance != -1)
write.table(P3F_C793_dep_data, "P3F_C793_dep_data.bed",col.names = F, row.names = F, quote = F, sep="\t")

P3F_C793_dependent_genelist<-unique(P3F_C793_dep_data$V14)


###############################################################################







