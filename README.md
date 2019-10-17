# DRIP-Seq Analysis - Raw data processing 
The dripseq_pipeline shell script takes fastq files from samples -/+OHT treatment.

After pre-processing, alignment and sorting, it uses a bed file to calculate read coverage around the cut asisi recognition sites.

This read coverage is normalised to the total read number to give reads per million for each nucleotide specified by the bed file.

Log2 fold change of +OHT(cut)/-OHT(uncut) can also be calculated to show the "damage-induced coverage".

The read coverage files can then be passed to the diva_plotting_pipe.R Rscript that generates a variety of plots.

# Plotting with diva_plotting_pipe.R

This script takes two positional arguments:
1. The directory in which the coverage files are saved.
2. The name of the coverage file.

All plots generated will be saved in the directory of the coverage file.

The following plots are generated:
1. Metagenes are plotted using a 10kb region centred on the break sites.
2. Boxplots of the same data as the metagenes are plotted by averaging the read coverage in a 1kb region centred on the break sites.
3. A correlation plot of coverage against transcriptional activity.
4. A heatmap of a 40kb region centred on the break sites and ordered from highest transcriptional activity to lowest. 

In addition, csv files are writted continaing the data used to generate the boxplots, the associated statistical tests and the correlation tests.

The Rscript can be edited to alter these parameters
