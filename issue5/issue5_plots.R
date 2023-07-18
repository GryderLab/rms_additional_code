# Yaw Asante | gryderlab | yxa181@case.edu | July 18th, 2023
# =============================================================================
# The following script uses RH4 HiChIP contacts and contextual p300, p3f, and cpg island
## data to show how contacts change at regions, segmented by whether they are overlapped
## by these features. This data shows how contacts mainly go down outside of CpG islands,
## but at regions bound by both p300 and P3F.


library(ggplot2)
library(cowplot)
library(tidyverse)

# using distance between Pol2 site and cpg island, here 0 to get overlap
DIST <- 0

# INPUT DATA
# ==============
# loads in a BED file of unified RH4_Pol2 regions between RH4_DMSO_6h_Pol2_HiChIP and
## RH4_dCBP1_6h_Pol2_HiChIP with L2FC of BAM coverage at each region, intersected with
## P3F binding sites, p300 binding sites, and cpg_islands
rh4_pol2_trip_ov <- read.delim("RH4_Pol2_P3F_p300_cpg_dist.l2fc.bed", header=F)

# DATA PROCESSING
# ==============
rh4_pol2_trip_ov <- rh4_pol2_trip_ov |> select(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V15)
colnames(rh4_pol2_trip_ov) <- c("chr","start","end","label","DMSO","dCBP1","diff","L2FC","wP3F","wP300", "distCPG")


# assigning labels to regions based on types of binding at each region
P3F_only <- which(rh4_pol2_trip_ov$wP3F > 0 & rh4_pol2_trip_ov$wP300 < 1)
P3F_p300 <- which(rh4_pol2_trip_ov$wP3F > 0 & rh4_pol2_trip_ov$wP300 > 0)
neither_P3F_p300 <- which(rh4_pol2_trip_ov$wP3F == 0 & rh4_pol2_trip_ov$wP300 == 0)
p300_only <- which(rh4_pol2_trip_ov$wP3F == 0 & rh4_pol2_trip_ov$wP300 > 0)
NEAR_CPG <- which(rh4_pol2_trip_ov$distCPG <= DIST)

rh4_pol2_trip_ov$group_type <- rep("neither_P3F_p300", dim(rh4_pol2_trip_ov)[1])
rh4_pol2_trip_ov$group_type[P3F_only] <- "Bound by P3F only"
rh4_pol2_trip_ov$group_type[p300_only] <- "Bound by p300 only"
rh4_pol2_trip_ov$group_type[P3F_p300] <- "Bound by p300 and P3F"
rh4_pol2_trip_ov$group_type[neither_P3F_p300] <- "Not bound by p300 or P3F"

rh4_pol2_trip_ov$CPG_SPLIT <- "No CpG Island Overlap"
rh4_pol2_trip_ov$CPG_SPLIT[NEAR_CPG] <- "Overlaps CpG Island"

# two types of data: type of p300 or P3F binding site and whether it overlaps a CpG island
rh4_pol2_trip_ov$group_type <- as.factor(rh4_pol2_trip_ov$group_type)
rh4_pol2_trip_ov$CPG_SPLIT <- as.factor(rh4_pol2_trip_ov$CPG_SPLIT)

# for clarity, we omit regions bound by P3F only (small group)
rh4_pol2_trip_ov <- rh4_pol2_trip_ov[rh4_pol2_trip_ov$group_type != "Bound by P3F only",]

# OUTPUT OF FINAL PLOT
# ==============
# writes the corresponding plot to a PDF for later polishing
pdf("RH4_Pol2_Contact_signal_by_P3F_p300_CpG.pdf", width=10, height=8)

ggplot(rh4_pol2_trip_ov, aes(y=L2FC, x=group_type)) + 
  geom_jitter(alpha=0.05, aes(color=group_type)) +
  geom_boxplot(width=0.5, outlier.shape = NA) + theme_minimal_grid() + 
  ggtitle("RH4 Pol2 Contact Signal") +
  xlab("Types of Regions by Occupancy") +
  ylab("L2FC of RNA Pol II Contacts") + facet_wrap(~ CPG_SPLIT)

graphics.off()

