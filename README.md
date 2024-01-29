# Additional Code for "PAX3-FOXO1 uses its activation domain to recruit CBP/P300 and shape RNA Pol2 cluster distribution"
Code provided by: Yaw Asante and Bhava Udhayakumar
Questions and data request can be made to: yaw.asante@case.edu

This repo contains the code and tools used to create the figures and datasets in Asante, Benischke and Osman et al, 2023. Any omissions or uncertainty can be reported to yaw.asante@case.edu


## Figures made, by Tool:

### ChIP-seq Track Visualization with IGV
Source: External, https://software.broadinstitute.org/software/igv/ 

Figure 3F, Figure 3H, Figure 5M

### Principal Component Analysis of RNA Seq
Source: In-House, in <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEA_Heat_PCA.R"> the RNA-seq pipeline in VisualizeRNAseq</a>

Figure 3A

### Ranked Gene Set Enrichment Analysis of Gene Expression
Source: In-House, <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist on GitHub</a>

Figure 3C, Figure 4H, Figure 5H

### Analysis of HiChIP Data
Source: In-House, using <a href="https://github.com/GryderLab/peaks3d">peaks3d</a> and aquatools.

Figures 6 (all panels), Figure 7 (all panels) and Figure 8 (all panels)



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

Panel F pairs profile plots of read-pileup with corresponding heatmaps and were made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/plotTSSheat">plotBEDheat</a> from the <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a> (which employs deeptools). 
Acetylation Levels and Chromatin Accessibility at p300 binding sites in RH4 cells.
<p align="center"><img width="472" alt="Metagene Plots for p300 and P3F binding" src="https://github.com/GryderLab/rms_additional_code/assets/135348829/dacf2277-1df0-4002-b61d-a829cf98c524">
</p>

Panel G is a principal component analysis of RNA-seq data, made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEA_Heat_PCA.R"> the RNA-seq pipeline </a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.

Panel H is a rank plot of gene set enrichment made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist script</a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.


### Figure 5
Panel H consists of rank plots of gene set enrichment made with <a href="https://github.com/GryderArt/VisualizeRNAseq/blob/master/RNAseq_Pipeline/buildTPM_Matrix_GSEAranklist_Heatmaps.R"> GSEA ranklist script</a> in <a href="https://github.com/GryderArt/VisualizeRNAseq/tree/master"> VisualizeRNASeq</a>.

Panel J is a bar plot made with GraphPad Prism.

Panel K is a line plot made with GraphPad Prism.

Panel M consists of ChIP-seq data tracks plotted with <a href="https://software.broadinstitute.org/software/igv/">Integrated Genome Viewer (IGV) </a>.

Panel N consists of profile plots of read-pileup with corresponding heatmaps and were made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/plotTSSheat">plotBEDheat</a> from the <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a> (which employs deeptools).


### Figure 6
Panel A (right) is output from the Genic Rank of Active Clustered Elements (GRACE) tool provided at <a href="https://github.com/GryderLab/peaks3d">peaks3d</a>. 

Panel B and Panel D were made in R with output from GRACE.

Panel E consists of boxplots made using RNA-seq and HiChIP data via code novel to the resumbission, able to be provided upon reasonable request.
Changes in RH4 gene expression due to dCBP1, split across features associated with genes near and distal from p300 enhancers. 
<p align="center"><img width="472" alt="Changes in expression due to dCBP1 with p300 + P3F and without P3F" src="https://github.com/GryderLab/rms_additional_code/assets/135348829/81d183e7-0452-4091-963a-44f36691f511">
</p> 

Panel F (left) are heatmaps and profile plots made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/plotTSSheat">plotBEDheat</a> from the <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a>.

Panel F (right) is an Aggregate Peak Analysis plot made with <a href="https://github.com/GryderLab/peaks3d">peaks3d</a>. 

Panel G consists of moving average plots made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/bedCovComp">BedCovCompare</a> from <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe repo</a>.

Panel H consists of rankplots of HiChIP data using BedCovCompare (as previously described here), created with code that can be provided upon reasonable request.
The impact of P300/CBP degraders on RH5 H3K27ac, segregated by enhancer and promoter regions (compared to impact on RH4 Pol2).
<p align="center"><img width="472" alt="BedCovCompare Output for Pol2 and H3K27ac HiChIP" src="https://github.com/GryderLab/rms_additional_code/assets/135348829/7126e44a-bafa-4f94-b832-4922e852269e">
</p>


### Figure 7
Panel A is a boxplot of HiChIP signal made with <a href="https://github.com/GryderLab/ChIPseqPipe/tree/master/bedCovComp">BedCovCompare</a> from <a href="https://github.com/GryderLab/ChIPseqPipe">ChIPseqPipe</a>.

Panel B was made using motif information from HOMER (version 4.9.1, http://homer.ucsd.edu/homer/index.html).

Panel C is a donut plot made with data from the HiChIP analysis obtained via <a href="https://github.com/GryderLab/peaks3d">peaks3d</a>.

Panel D is a boxplot made using same HiChIP data from Panel C overlapped with p300 and P3F binding sites from ChIP-seq as described, with code which can be provided upon reasonable request.
3D genome looping at Pol2 sites strongly bound by P300 with and without P3F.
<p align="center"><img width="472" alt="Looping to sites with and without P3F" src="https://github.com/GryderLab/rms_additional_code/assets/135348829/7ebd86d5-3b27-4b96-8351-b694422ef3e7">
</p>

Panel E is a histogram made using the same HiChIP data from Panel C.

### Figure 8
Panel B made using HiChIP data viewed via <a href="https://github.com/aidenlab/Juicebox">Juicebox</a>, <a href="https://software.broadinstitute.org/software/igv/">Integrated Genome Viewer (IGV) </a> and analyzed with plotAPA and runpeaks3d (generation of loops) from <a href="https://github.com/GryderLab/peaks3d">peaks3d</a>.

Panel D consists of Aggregate Peak Analysis (APA) plots of RH4 HiChIP data (with dTAG47 KO of P3F) made with plotAPA from <a href="https://github.com/GryderLab/peaks3d">peaks3d</a>. APA plots made with the parameters described in the figure, but code can be provided upon reasonable request.
