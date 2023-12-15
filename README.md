# Additional Code for "PAX3-FOXO1 uses its activation domain to recruit CBP/P300 and shape RNA Pol2 cluster distribution"
by Yaw Asante and Bhava Udhayakumar <br>
(Updated: 12-15-2023 -> Usability and readability updates coming later this week, thanks for your early interest! Immediate questions can be sent to yaw.asante@case.edu) <br>

This repo contains the code and tools used to create the figures and datasets in Asante, Benischke and Osman et al, 2023. Other tools used indicated as noted at the bottom.

## Figures made, by Tool:

### ChIP-seq Track Visualization with IGV
Source: External, https://software.broadinstitute.org/software/igv/ 
Figure 3F, Figure 3H, Figure 5M

### Principal Component Analysis of RNA Seq
Source: In-House, in 
Figure 3A

### Ranked Gene Set Enrichment Analysis of Gene Expression
Source: In-House, <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist on GitHub</a>
Figure 3C, Figure 4H, Figure 5H

### Analysis of HiChIP Data
Source: In-House
Figures 6, 7 and 8

## Tools used by Figure:

### Figure 3
Panel A is a principal component analysis of RNA-seq derived TPM data, made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEA_Heat_PCA.R"> the RNA-seq pipeline </a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.
Panel B is a heatmap of expression made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEA_Heat_PCA.R"> the RNA-seq pipeline </a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.
Panel C is a rank plot of gene set enrichment made with the <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist script</a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.
Panel D is a scatterplot made with GraphPad Prism.
Panel E is a pie chart made with GraphPad Prism.
Panel F and Panel H consist of ChIP-seq data tracks plotted with the <a href="https://software.broadinstitute.org/software/igv/">Integrated Genome Viewer (IGV) </a>.
Panel G is a bar plot made with GraphPad Prism. 
Panel I pairs a profile plot and a heatmap and was made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/plotTSSheat">plotBEDheat</a> from the <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a> (which employs deeptools).

### Figure 4
Panel A contains a volcano plot made with GraphPad Prism.
Panel F pairs profile plots of read-pileup with corresponding heatmaps and were made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/plotTSSheat">plotBEDheat</a> from the <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a> (which employs deeptools).. 
Panel G is a principal component analysis of RNA-seq data, made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEA_Heat_PCA.R"> the RNA-seq pipeline </a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>..
Panel H is a rank plot of gene set enrichment made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist script</a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.

### Figure 5
Panel H consists of rank plots of gene set enrichment made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist script</a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.
Panel J is a bar plot made with GraphPad Prism.
Panel K is a line plot made with GraphPad Prism.
Panel M consist of ChIP-seq data tracks plotted with <a href="https://software.broadinstitute.org/software/igv/">Integrated Genome Viewer (IGV) </a>.

### Figure 6
Panel A (right) is output from the Genic Rank of Active Clustered Elements (GRACE) tool provided at 
Panel B and Panel D were made with output from GRACE.
Panel E consists of boxplots made using code novel to the resumbission, provided in here in subfolder 
Panel F (left) are heatmaps and profile plots made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/plotTSSheat">plotBEDheat</a> from the <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a>.
Panel F (right) is an Aggregate Peak Analysis plot made with <a href=""> </a> 
Panel G consists of moving average plots made with <a href=""> </a> from <a href=""> </a>.
Panel H consists of rankplots of HiChIP data using data novel to the resubmission, provided here in subfolder 6H.

### Figure 7
Panel A is a boxplot of HiChIP signal made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/bedCovComp">BedCovCompare</a> from <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe</a>.
Panel B was made using motif information from HOMER (version 4.9.1, http://homer.ucsd.edu/homer/index.html).
Panel D is a boxplot made using HiChIP data with code novel to the resubmission, provided here in subfolder 7D.

### Figure 8
Panel D consists of Aggregate Peak Analysis (APA) plots made with plotAPA from <a href="https://github.com/GryderLab/peaks3d">peaks3d</a>. Data used was novel to the resubmisison.

### Issue 1
Activity from directly bound P3F target genes (within 500kb of P3F binding) vs indirect targets. 
(Code provided, data included as part of Supplementary Table 2, model plot included below)
<p align="center">![image](https://github.com/GryderLab/rms_additional_code/assets/135348829/d400b390-c94c-41e0-ac9c-1685462a612d)
</p>
(Genes in direct: PIPOX, FRG2B, DBX1, PERCC1, CDH4, KIF3B, FGF8, FOXF1, FENDRR, EFNB2, FGFR2)

### Issue 2
Separation into C793 mutation dependent and partially dependent gene sets combined with an analysis of the contextual nature of the P3F binding site motifs.
(Output from HOMER provided, no meaningful differences)
<p align="center">  
</p>

### Issue 3
Acetylation Levels and Chromatin Accessibility at p300 binding sites in RH4 cells.
<p align="center">![image](https://github.com/GryderLab/rms_additional_code/assets/135348829/dacf2277-1df0-4002-b61d-a829cf98c524)
</p>

### Issue 4
The impact of P300/CBP degraders on RH5 H3K27ac, segregated by enhancer and promoter regions (compared to impact on RH4 Pol2).
<p align="center">![image](https://github.com/GryderLab/rms_additional_code/assets/135348829/7126e44a-bafa-4f94-b832-4922e852269e)
</p>
<p align="center"><img width="472" alt="MA plots for H3K27ac at TSS proximal and distal regions" src="https://github.com/GryderLab/rms_additional_code/assets/135348829/e999d65c-282e-458a-beef-25254d3dd454">
</p>

### Issue 5
3D genome looping at Pol2 sites strongly bound by P300 with and without P3F.
<p align="center">![image](https://github.com/GryderLab/rms_additional_code/assets/135348829/7ebd86d5-3b27-4b96-8351-b694422ef3e7)
</p>

### Issue 6
Changes in RH4 gene expression due to dCBP1, split across genes near and distal from p300 enhancers.
<p align="center">![image](https://github.com/GryderLab/rms_additional_code/assets/135348829/81d183e7-0452-4091-963a-44f36691f511)

</p>
