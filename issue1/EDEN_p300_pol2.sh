#!/bin/bash
#SBATCH -n 1
#SBATCH -p smp 
#SBATCH --mem=4g 
#SBATCH --cpus-per-task=8 
#SBATCG --time=1:00:00

##################################################
# Eden pipeline SG 2021 -11-21     
##################################################

############################################
# Declaring variables
HPC='/mnt/rstor/SOM_GENE_BEG33'

HOME_DNA=${HPC}/ChIP_seq/hg38/DATA
HOME_RNA=${HPC}/RNA_seq/DATA

############################################
# loading modules
module load miniconda3/4.9.2
source activate "/home/sxg1131/.conda/envs/py368"

# -t Cufflink expression: 'EXP' , Khanlab TPM: 'TPM', Cuffdiff: 'DIFF', Matrix file: 'MATRIX'
# -f  <float>   TPM cutoff (default: $TPM_cutoff)
# -b  <string>  Input BED file 
# -e  <string>  expression file (UNIX format!)
# -d  <string>  TAD file
# -o  <string>  Output directory (default: /mnt/rstor/SOM_GENE_BEG33/sxg1131/testEden_20211121)
# -x  <string>  Output file prefix
# -n  <integer> The distance to nearest gene (default: $nearest_gene_loci_cutoff)

#run EDEN for RH4 with MAX TPM to find genes
perl EDEN.pl -b RH4_P300_enh_w_P3F.bed -e RH4_p300_Pol2.coding.TPM.matrix.txt -d TAD_goldmine.hg38.bed -o /mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2_real/P3F_motifs -x P3F_w_Pol2_genelist -t TPM -c -f 1 -n 1000000 "4,5,6,7,8"

perl EDEN.pl -b RH4_P300_enh_no_P3F.bed -e RH4_p300_Pol2.coding.TPM.matrix.txt -d TAD_goldmine.hg38.bed -o /mnt/rstor/SOM_GENE_BEG33/ChIP_seq/hg38/projects/p300_Pol2_real/P3F_motifs -x P3F_no_Pol2_genelist -t TPM -c -f 1 -n 1000000 "4,5,6,7,8"


