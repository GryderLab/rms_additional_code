# rms_additional_code
This repo contains the code used to create the figures and datasets requested by the reviewers of the Osman et al 2023 paper, submitted to Nature Communications.

1. Due to the timing of the experiment it is difficult to separate direct from indirect targets. If available, the authors should look at directly bound P3F target genes (data from Figure 4F). 

We developed box plots which showed, for these genes, which ones have P3F bound in a nearby enhancer.

2. The authors only use an FDR cutoff in their analysis, but a cutoff taking into consideration the magnitude of the effect would be important. Moreover, separation into most sensitive and least sensitive genes combined with an analysis of the contextual nature of the P3F binding sites may reveal whether certain genes are more or less affected by the C793 mutation and possibly why.

We demonstrate this by comparing the motifs at C793 dependent and partially dependent genes (HOMER).

3. The authors demonstrate that p300 binding is enriched at sites occupied by P3F, how does this compare with acetylation levels and chromatin accessibility at these sites?


4. The impact of P300i and P300-PROTACs on histone acetylation H3K18ac or H3K27ac should be assessed and segregated by enhancer and promoter regions. 

To demonstrate this by using our RH5 H3K27ac data, and plotting the read signal at the sites most altered by dCBP1 treatment, split by whether they lie proximal (<= 5kb) or distal (> 5kb) from gene promoters.
(BedCovCompare, moving average plots.png)

5. What happens to 3D genome looping at Pol2 sites strongly bound by P300 but that are not occupied by P3F. Are they affected in a similar manner? 

6. The models in Figure 7D-E are based on the results presented previously (Figure 6F, G, H and I and Figure 7 A-B), it is not clear to this reviewer that CBP/P300 inhibition-induced migration of the RNApol 2 cluster from enhancers to promoter-enriched CpGs preferentially in P3F target genes, beyond the examples shown for SOX8 and MYOD1 (Figure 7A-B). Controls that support the proposed hypothesis are required.
7. CBP/P300 are central/core regulators of gene activation in cells. How expected are the results of the second part? Can the authors put more emphasis what’s going on outside PAX3-FOXO1 target regions and discuss more this.
