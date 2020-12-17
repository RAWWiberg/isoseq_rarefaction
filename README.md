# cDNA_Cupcake pipeline: 
### isoseq3 post-processing and rarefaction analysis

An obvious question to ask oneself is "how much sequencing do I need to do"? With IsoSeq being a relatively new technology, and quite expensive it would be good to know how many more transcripts/genes one can recover with additional sequencing. One procedure for figuring this out is a rarefaction analysis where one downsamples from the number of reads one has produce to see how many are recovered, as the number of reads increases. One such solution is implemented in the popular cDNA Cupcake ToFU pipeline: [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake).

However, in this pipeline it's not clear that this rarefaction analysis is entirely satisfying. Because the "standard" used for the analysis is all the discovered transcripts in the full sample, it is really only measuring how many additional transcripts will be discovered with more sequencing effort. Clearly there are not an infinite number of transcripts available to sample within a cell. So perhaps a more sensible question is "how much sequencing effort do we need to get a good representation of full complement of transcripts within the sample?"  

To answer this question a different standard is needed, a standard that is, at least to some degree, a measure of the total complement of transcripts/genes.

We performed a very similar rarefaction analysis as what is done in cDNA_Cupcake, but instead of measuring how many of the set of discovered transcripts are present in each subsample, we measured the percent (%) of complete BUSCO (duplicated and single-copy) genes that are present in each subsample.

Full details can be followed in the R script that you can find in the `scripts` folder of this repository.

Briefly, for this procedure I use the clusters produced by the isoseq3 pipeline along with the `*cluster_report.csv` file which is produced bu the `isoseq3 cluster` step.
I also use the file called `full_table.tsv` which is produced by running BUSCO on the final clustered transcripts (found in the file `*.clustered.hq.fasta.gz` after running the `isoseq3 cluster` step.

The `cluster_report.csv` file gives the ID of each identified cluster as well as which CCS reads support that ID. I can use this to downsample reads, set a threshhold for the number of reads I need to support an ID, and ask how many of the reported IDs I recover when I downsample to X reads. In the file `full_table.tsv` from the BUSCO output each cluster ID is mapped to a BUSCO gene ID. I can use this file to then ask of the recovered IDs at X downsampled reads, how many BUSCO genes would I be recovering. This becomes my Y variable for each increment of X (the number of CCS reads). For each increment of X I do the downsampling Z times to get a distribution of values.

![fig1][fig1]

The figure above shows the mean (solid black line), SE (black dashed lines), and max (dashed red line) percentage of BUSCO genes recovered as a function of the number of CCS reads sampled. For the above figure X takes a value between 0 and 199000 (in increments of 1000) and I do 100 re-sampling rounds for each increment of X.

Also shown are some benchmark values of the percentage of BUSCO genes recovered from 1) a whole-genome assembly of the same organism from which our Iso-Seq data come, 2) a genome-guided transcriptome assembly based on the genome assembly and Illumina RNA-seq data, and 3) the isoseq3 clusters for reference.

We can zoom in on the data for samples with >200,000 reads to see that the trend is almost linear at this point (though not quite).

![fig2][fig2]

If we fit a regression line to that data, and extrapolate out to how many reads we need to cross 100% BUSCO genes, we obtain the following figure.

![fig3][fig3]

Where the grey dashed line charts the extrapolations from the regression model. Only at ~1.5-1.8 million CCS reads does it cross the threshhold for the genome assembly.

We can test this intuition if we produce more Iso-Seq data. Which is what I did here. I obtained the above results with a sample that was run on a single 1M SMRT cell on the Sequel system. This produced 200,662 CCS reads represented the final set of high quality transcripts. We subsequently ran the same sample on an 8M SMRT cell on the Sequel II system to produce an total of 923,429 CCS (i.e. an additional 772,767 CCS reads). So I now run the same procedure as above and ask how many BUSCO genes do I recover using the full set of ~920000 CCS reads. 

From the figure above this should put us at ~50% completeness. 

![fig4][fig4]


From the figure above we see that we perform a bit worse than expected (red point).
In reality, we perform a bit worse than expected since a linear regression model is not 100% appropriate in this case, there is still some curvature to the trend. But...we do still gain more BUSCO genes which suggests there is still some diversity in the IsoSeq library that we did not initially capture, which warrants further sequencing.

However, we are clearly doing better than we would expect had we fit a curved line to the initial rarefaction data.

![fig5][fig5]

In the figure above I have performed the same rarefaction analysis with all data combined (blue solid line) and compare it to the same analysis with just the first sample (red solid line).
If we extended the curved red line further to the right, we would not expect to see the improvement we did. 
Thus, it's difficult to make good predictions from a subset of a library about the number of transcripts that are in fact present.

### Conclusions

In our case one SMRT cell was clearly not enough to get a good representation of the genes. 

Based on the extrapolations we would need >2-3 million CCS reads to achieve the same results as the genome assembly threshhold


#### Data statement:
All of the data and scripts required to check these observations are included in this repository.


[fig1]:/figures/sp1_plot1.png
[fig2]:/figures/sp1_plot2.png
[fig3]:/figures/sp1_plot1_plus1.png
[fig4]:/figures/sp1_plot1_plus2.png
[fig5]:/figures/sp1_plot1_plus3.png


