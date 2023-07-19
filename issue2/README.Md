Author: Yaw Asante | yxa181@case.edu | Gryder Lab
Date Created: July 18th, 2023

Comment from Reviewer:
For genes which have P3F bound in a nearby enhancer, compare these P3F peaks to P3F peaks that are not near a C793 sensitive gene.  Are the motifs at that peak different?
> Response: Using HOMER, we identified and plotted the binding site motifs of those peaks associated with dependent genes and those which were not.

Creating the Input Data File(s):
To identify the two sets of P3F peaks (those associated with failed P3F_C793S rescue and those associated with partial rescue), we used bedtools to intersect the transcription start sites of the genes from either group (mutation dependent or mut_dep and mutation partially dependent or par_dep) with P3F binding sites identified from ChIP_seq data from GSE116344 (MACS threshold p-14).


We then ran findMotifsGenome.pl (from GryderChIPseq_scripts, implemented with homer version 4.11) on these BED files to get the motifs associated with either group.

Creating the Primary Output Figure:
All of the code required for making the plot is in: issue2_plots.R

Final Figure:
The final figure was created by altering label placement from the primary output figure in Adobe Illustrator.
