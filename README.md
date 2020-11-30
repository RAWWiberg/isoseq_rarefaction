# cDNA_Cupcake pipeline: 
### isoseq3 post-processing and rarefaction analysis

[cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake)

However, it's not clear what this rarefaction analysis is aiming to accomplish. Because the "standard" used for the analysis is all the discovered transcripts, it is really only measuring how many additional transcripts will be discovered with more sequencing effort. Clearly there are not an infinite number of transcripts available to sample within a cell. So perhaps a more sensible question is "how much sequencing effort do we need to get a good representation of full complement of transcripts within the sample?"  

To answer this question a different standard is needed, a standard that is, at least to some degree, a measure of the total complement of transcripts.

We performed a very similar rarefaction analysis as what is done in cDNA_Cupcake, but instead of measuring how many of the set of discovered transcripts are present in each subsample, we measured the percent (%) of complete BUSCO (duplicated and single-copy) genes that are present in each subsample.

Full details can be followed in the R script that you can find in the `scripts` folder of this repository.

Briefly, I use the clusters produced by the isoseq3 pipeline along with the `*cluster_report.csv` file which is produced bu the `isoseq3 polish` step


![fig1][fig1]

The figure above shows the mean (solid black line), SE (black dashed lines), and max (dashed red line) percentage of BUSCO genes recovered as a function of the number of CCS reads sampled.
Also shown are some benchmark values of the percentage of BUSCO genes recovered from a whole-genome assembly of the same organism from which our Iso-Seq data come, a genome-guided transcriptome assembly based on the genome assembly and Illumina RNA-seq data, and the isoseq3 clusters for reference.

We can zoom in on the data for samples with >200,000 reads to see that the trend is almost linear at this point (though not quite).

![fig2][fig2]

If we fit a regression line to that data, and extrapolate out to 4 million CCS reads. We obtain the following figure.

![fig3][fig3]

Where the grey dashed line charts the extrapolations from the regression model. Only at ~4 million CCS reads does it cross the threshhold for the genome assembly.
In reality, it will probably cross much later since a linear regression model is not 100% appropriate in this case, there is still some curvature to the trend.

### Conclusions

In our case one SMRT cell was clearly not enough to get a good representation of the genes. 

Based on the extrapolations we would need ~4 million CCS reads to achieve the same results as the genome assembly threshhold. This corresponds to ~16 SMRT cells on the Sequel and 2-3 on the Sequel II.



#### Data statement:
All of the data and scripts required to check these observations are included in this repository.


[fig1]:/figures/Maccli_plot1.png
[fig2]:/figures/Maccli_plot2.png
[fig3]:/figures/Maccli_plot1_plus.png



