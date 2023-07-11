# Motif L2FC

mut_dep <- read.delim("mut_dep_knownResults.txt", sep="\t")
mut_indep <- read.delim("indep_knownResults.txt", sep="\t")


mut_dep <- mut_dep[, c(1,7,9)]
mut_indep <- mut_indep[, c(1,7,9)]

colnames(mut_dep)[c(2,3)] <- c("dep_occ", "dep_bkg")
colnames(mut_indep)[c(2,3)]<- c("indep_occ","indep_bkg")

mot_table <- left_join(mut_dep, mut_indep, by="Motif.Name")

mot_table[,2] <- gsub("%","", mot_table[,2])
mot_table[,3] <- gsub("%","", mot_table[,3])
mot_table[,4] <- gsub("%","", mot_table[,4])
mot_table[,5] <- gsub("%","", mot_table[,5])

mot_table[,2] <- as.double(mot_table[,2])
mot_table[,3] <- as.double(mot_table[,3])
mot_table[,4] <- as.double(mot_table[,4])
mot_table[,5] <- as.double(mot_table[,5])
mot_table$L2FC <- log((mot_table$dep_occ / mot_table$dep_bkg) / (mot_table$indep_occ/mot_table$indep_bkg))
#mot_table <- mot_table[order(mot_table$L2FC, decreasing=T),]
#mot_table <- mot_table[order(mot_table$L2FC, decreasing=T),]

mot_table$wL2FC <- mot_table$L2FC * mot_table$dep_occ
mot_table <- mot_table[order(mot_table$wL2FC, decreasing=T),]
mot_table$rank <- c(1:dim(mot_table)[1])

ggplot(mot_table, aes(x=rank, y=wL2FC)) + geom_point()


mot_top <- mot_table[mot_table$wL2FC > 0,]
ggplot(mot_top, aes(x=rank, y=wL2FC)) + geom_point()

