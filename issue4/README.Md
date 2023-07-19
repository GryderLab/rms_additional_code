Author: Yaw Asante | yxa181@case.edu | Gryder Lab
Date Created: July 18th, 2023

## Comment from Reviewer:
The impact of P300i and P300-PROTACs on histone acetylation H3K18ac or H3K27ac should be assessed and segregated by enhancer and promoter regions. 
> Response: Using the phenotypically similar RH5 cell-line, we assessed the impact of p300 degradation on HiChIP contacts in H3K27ac data by comparing the read coverage at the union of peaks between the control (DMSO) and dCBP1 cases, split between regions proximal (<= 5kb) and distal (> 5kb) from gene promoters.

## Creating the Input Data File(s)
The process of running this data was as described in https://github.com/GryderLab/ChIPseqPipe/tree/master/bedCovComp for RH4 RNA Pol2 HiChIP and RH5 H3K27ac HiChIP, in both cases comparing the DMSO- and dCBP1-treated cells at 6 hours.

## Creating the Primary Output Figure
The invoke used to create the data shown here for RH4 is:
```
sbatch -p smp --ntasks=2 --mem=32g --cpus-per-task=16 -J bCComp --export=sample1=RH4_DMSO_6h_Pol2_HiChIP_081721_NovoG,sample2=RH4_dCBP1_6h_Pol2_HiChIP_081721_NovoG,name1=RH4_DMSO_6h_Pol2,name2=RH4_dCBP1_6h_Pol2,pValue=p-4,overlapname=gene_TSS5k,overlapbed=gene_tss5k.clean.bed runBedCovCompareSpiked.sh
```
The invoke used to create the data shown here for RH5 is:
```
sbatch -p smp --ntasks=2 --mem=32g --cpus-per-task=16 -J bCComp --export=sample1=RH5_DMSO_H3K27ac_070622_NovoG,sample2=RH5_dCBP1_H3K27ac_070622_NovoG,name1=RH5_DMSO_H3K27ac,name2=RH5_dCBP1_H3K27ac,pValue=p-4,overlapname=gene_TSS5k,overlapbed=gene_tss5k.clean.bed runBedCovCompareSpiked.sh
```
> Note that in both cases, the script run is "runBedCovCompareSpiked.sh", which runs the comparison informed by AQuA values (ratios between the hg38 reads and the orthogonal mm10 reads sequenced at the same time).


## Final Figure
The final figure was created by adjusting the axes, removing extraneous labels, and placing side-by-side the primary output figures, alongside accompanying cartoons of what appears to be taking place (created with Adobe Illustrator).
