Author: Yaw Asante | yxa181@case.edu | Gryder Lab
Date Created: July 18th, 2023

Comment from Reviewer:
What happens to 3D genome looping at Pol2 sites strongly bound by P300 but that are not occupied by P3F?
> Response: We assessed this by using the L2FC of read coverage between our RH4 RNA Polymerase-II targeted HiChIP at the union of peaks between control (DMSO) and drug-treated (dCBP1) cases, splitting the regions by whether they overlap p300 only, P3F only, or p300 and P3F binding sites. This is further split by CpG island occupancy.

Creating the Input Data File(s)
The initial BED file was created using BedCovCompare (see: https://github.com/GryderLab/ChIPseqPipe/tree/master/bedCovComp) on RH4 Pol2 HiChIP data comparing BAM coverage of the DMSO and dCBP1 treated cells after 6 hours. For context, I've included the primary output file of BedCovCompare used here.

sbatch -p smp --ntasks=2 --mem=32g --cpus-per-task=16 -J bCComp --export=sample1=RH4_DMSO_6h_Pol2_HiChIP_081721_NovoG,sample2=RH4_dCBP1_6h_Pol2_HiChIP_081721_NovoG,name1=RH4_DMSO_6h_Pol2,name2=RH4_dCBP1_6h_Pol2,pValue=p-4,overlapname=gene_TSS5k,overlapbed=gene_tss5k.clean.bed runBedCovCompareSpiked.sh
> Gives you RH4_DMSO_6h_Pol2_RH4_dCBP1_6h_Pol2_gene_TSS5k.overlap.AB.delta.l2fc.bed

bedtools intersect -a RH4_DMSO_6h_Pol2_RH4_dCBP1_6h_Pol2_gene_TSS5k.overlap.AB.delta.l2fc.bed -b RH4_D6_PAX3FOXO1_GSE116344.nobl.p-14.bed -c > RH4_Pol2_P3F_overlap.l2fc.bed #MACS threshold 0.00000000000001
bedtools intersect -a RH4_Pol2_P3F_overlap.l2fc.bed -b RH4_p300SC_GSE83726.nobl.bed -c > RH4_Pol2_P3F_p300_overlap.l2fc.bed #MACS threshold 0.00001
bedtools intersect -a RH4_Pol2_P3F_p300_overlap.l2fc.bed -b cpgisland.hg38.bed > RH4_Pol2_P3F_p300_cpg_dist.l2fc.bed

Creating the Primary Output Figure
All of the code required for making the plot is in: issue5_plots.R

Final Figure
The final figure was created by altering label placement from the primary output figure in Adobe Illustrator.
