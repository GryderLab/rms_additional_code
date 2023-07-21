getwd()
setwd("/Volumes/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2_real/P3F_motifs")

library(plyr)

# setwd("/Volumes/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2/P3F_motifs")
##lists of genes that are rescued or not by C793S (from the paper) 
###two gene lists

C793S_rescued<-read.delim("C793S_dependent_genes.txt", header = F, sep = " ")
C793S_rescued<-C793S_rescued$V1
C793S_partial_rescue<-read.delim("C793S_partial_rescue_independent_genes.txt", header = F, sep = " ")
C793S_partial_rescue<-C793S_partial_rescue$V1

#get the gene info from TSS bed file
Sample_Name = "P3F"
Folder_Name = "P3F_motifs"
GENE_bed<-read.table("gene_tss5k.clean.bed")
colnames(GENE_bed)<-c("Chrom", "Start", "End", "Gene", "Score", "Strand")

##filtering for the C793S rescued or partial rescued genes from the TSS file
###C793_rescued
C793S_rescued.genetable <- GENE_bed[GENE_bed$Gene %in% C793S_rescued, ]
missing_genes_1 <- C793S_rescued[!(C793S_rescued %in% C793S_rescued.genetable$Gene)]

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
#write.table(C793S_rescued.genetable, "TSS_C793S_rescued_genelist.bed",col.names = F, row.names = F, quote = F, sep="\t")

###C793_partial_rescue
C793S_partial_rescue.genetable<-GENE_bed[GENE_bed$Gene %in% C793S_partial_rescue, ]
missing_genes_2 <- C793S_partial_rescue[!(C793S_partial_rescue %in% C793S_partial_rescue.genetable$Gene)]

# # Identify genes repeated twice
# repeated_genes <- C793_low.dependent.genetable %>%
#   group_by(Gene) %>%
#   filter(n() == 2) %>%
#   distinct(Gene)
# repeated_genes$Gene
#ICOSLG gene repeated twice
#write.table(C793S_partial_rescue.genetable, "TSS_C793S_partial_rescue.genelist.bed",col.names = F, row.names = F, quote = F, sep="\t")

###############################################################################
#filtering for genes within 500kb from the bedtools closest output
##C793 partial rescue
RH4_C793S_partial_rescue<-read.table("P3F_near_C793_low.dependent.bed")

##rename bedtools closest col names
# colnames(P3F_C793dependent) <- c("chrom_A", "start_A", "end_A", "name_A", "score_A", "strand_A", "chrom_B", "start_B", "end_B", "name_B", "score_B", "strand_B", "distance")

req_distance=500000
RH4_C793S_partial_rescue$distance=RH4_C793S_partial_rescue[,ncol(RH4_C793S_partial_rescue)]

RH4_direct_C793S_partial_rescue <- subset(RH4_C793S_partial_rescue, RH4_C793S_partial_rescue$distance >= -(req_distance) & RH4_C793S_partial_rescue$distance <= req_distance)

RH4_indirect_C793S_partial_rescue <- subset(RH4_C793S_partial_rescue, RH4_C793S_partial_rescue$distance < -(req_distance) | RH4_C793S_partial_rescue$distance > req_distance)

RH4_direct_C793S_partial_rescue <- subset(RH4_direct_C793S_partial_rescue, RH4_direct_C793S_partial_rescue$distance != -1)

write.table(RH4_direct_C793S_partial_rescue, "P3F_C793S_direct_partial_rescue_data.bed",col.names = F, row.names = F, quote = F, sep="\t")

RH4_direct_C793S_partial_rescue_genelist<-unique(RH4_direct_C793S_partial_rescue$V14)
RH4_indirect_C793S_partial_rescue_genelist<-unique(RH4_indirect_C793S_partial_rescue$V14)

##C793 dependent
RH4_C793S_rescued<-read.table("P3F_near_C793_dependent.bed")

##rename bedtools closest col names
# colnames(P3F_C793dependent) <- c("chrom_A", "start_A", "end_A", "name_A", "score_A", "strand_A", "chrom_B", "start_B", "end_B", "name_B", "score_B", "strand_B", "distance")

RH4_C793S_rescued$distance=RH4_C793S_rescued[,ncol(RH4_C793S_rescued)]

RH4_direct_C793S_rescued <- subset(RH4_C793S_rescued, RH4_C793S_rescued$distance >= -(req_distance) & RH4_C793S_rescued$distance <= req_distance)
RH4_indirect_C793S_rescued <- subset(RH4_C793S_rescued, RH4_C793S_rescued$distance < -(req_distance) | RH4_C793S_rescued$distance > req_distance)

RH4_direct_C793S_rescued <- subset(RH4_direct_C793S_rescued, RH4_direct_C793S_rescued$distance != -1)
write.table(RH4_direct_C793S_rescued, "P3F_C793S_direct_rescued_data.bed",col.names = F, row.names = F, quote = F, sep="\t")

RH4_direct_C793S_rescued_genelist<-unique(RH4_direct_C793S_rescued$V14)
RH4_indirect_C793S_rescued_genelist<-unique(RH4_indirect_C793S_rescued$V14)
###############################################################################

##4 genelists
RH4_direct_C793S_partial_rescue_genelist
RH4_indirect_C793S_partial_rescue_genelist
RH4_direct_C793S_rescued_genelist
RH4_indirect_C793S_rescued_genelist

##pie chart
gene.distribution.groups <- matrix(c(299, 351, 56, 79), nrow = 2, ncol = 2, dimnames = list(c("Direct", "Indirect"),
                                      c("C793S rescued", "C793S partial rescue")))

labels <- c("Direct C793S rescued: 299", "Indirect C793S rescued: 351",
           "Direct C793S partial rescue: 56", "Indirect C793S partial rescue: 79")
percent_labels <- paste0(labels, " (", round(gene.distribution.groups/sum(gene.distribution.groups) * 100), "%)")
pie(gene.distribution.groups, main = "Distribution within the groups", labels = percent_labels)












