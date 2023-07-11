library(ggplot2)
library(cowplot)
library(tidyverse)

# using distance between Pol2 site and cpg island (less than DIST threshold)
DIST <- 0

rh4_pol2_trip_ov <- read.delim("RH4_Pol2_P3F_p300_cpg_dist.l2fc.bed", header=F)
rh4_pol2_trip_ov <- rh4_pol2_trip_ov |> select(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V15)
colnames(rh4_pol2_trip_ov) <- c("chr","start","end","label","DMSO","dCBP1","diff","L2FC","wP3F","wP300", "distCPG")

#rh4_pol2_trip_ov$P3F_only <- (rh4_pol2_trip_ov$wP3F > 0 & rh4_pol2_trip_ov$wP300 < 1)
#rh4_pol2_trip_ov$P3F_p300 <- (rh4_pol2_trip_ov$wP3F > 0 & rh4_pol2_trip_ov$wP300 > 0)


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

rh4_pol2_trip_ov$CPG_SPLIT <- "Distant from CpG"
rh4_pol2_trip_ov$CPG_SPLIT[NEAR_CPG] <- "Near CpG"
#rh4_pol2_trip_ov$CPG_SPLIT <- "No CpG overlap"
#rh4_pol2_trip_ov$CPG_SPLIT[NEAR_CPG] <- "CpG overlap"


rh4_pol2_trip_ov$group_type <- as.factor(rh4_pol2_trip_ov$group_type)
rh4_pol2_trip_ov$CPG_SPLIT <- as.factor(rh4_pol2_trip_ov$CPG_SPLIT)

# BK Requested corrections
rh4_pol2_trip_ov <- rh4_pol2_trip_ov[rh4_pol2_trip_ov$group_type != "Bound by P3F only",]

rh4_pol2_trip_ov$STATUS <- rep(FALSE, dim(rh4_pol2_trip_ov)[1])
rh4_pol2_trip_ov$STATUS <- rep(FALSE, dim(rh4_pol2_trip_ov)[1])

# for conversion to count by split tables
rh4_pol2_trip_ov_dn <- rh4_pol2_trip_ov[rh4_pol2_trip_ov$L2FC <= 0,]
rh4_pol2_trip_ov_up <- rh4_pol2_trip_ov[rh4_pol2_trip_ov$L2FC > 0,]

up_split <- table(rh4_pol2_trip_ov_up$group_type, rh4_pol2_trip_ov_up$CPG_SPLIT)

pdf("RH4_Pol2_Contact_signal.pdf", width=10, height=8)
ggplot(rh4_pol2_trip_ov, aes(y=L2FC, x=group_type)) + 
  #stat_ydensity(geom="violin", aes(fill=group_type)) + 
  geom_jitter(alpha=0.05, aes(color=group_type)) +
  geom_boxplot(width=0.5, outlier = NA) + theme_minimal_grid() + 
  #theme(text=element_text(size=16)) +
  ggtitle("RH4 Pol2 Contact Signal") +
  xlab("Types of Regions by Occupancy") +
  ylab("RNA Pol II Log2FC") + facet_wrap(~ CPG_SPLIT)

graphics.off()

rh4_ext_trip <- rh4_pol2_trip_ov[(rh4_pol2_trip_ov$L2FC) < -2,]
ggplot(rh4_ext_trip, aes(y=L2FC, x=group_type)) + 
  #stat_ydensity(geom="violin", aes(fill=group_type)) + 
  geom_point(alpha=0.3, aes(color=group_type)) +
  #geom_boxplot(width=0.1) +
  theme_minimal_grid() + 
  #theme(text=element_text(size=16)) +
  ggtitle("Regions in RH4 with Strong L2FC DN") +
  xlab("Types of Regions by Occupancy") +
  ylab("RNA Pol II Log2FC")



data <- data.frame(group=factor(c("Bound by P3F only (21)", "Bound by p300 only (16464)", "Bound by p300 and P3F (1886)",
                                  "Not bound by p300 or P3F (16747")),
                   value=c(21, 16464, 1886, 16747))
#ggplot(count_table, aes(x="", y=number, fill=case)) +
#  geom_bar(stat="identity", width=1) +
#  coord_polar("y", start=0)

data <- data %>% 
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

print(data)

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  
  geom_text(aes(y = ypos, label = group), color = "white", size=6) +
  scale_fill_brewer(palette="Set1")
