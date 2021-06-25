# DRIP-Seq Analysis - Raw data processing 

This pipeline is for use with drip_seq data generated at DNA double-strand breaks induced by the DIvA (Damaged Induced via AsiSI) cell system.

The dripseq_pe_pipeline shell script takes paired end fastq files/bcl files from samples -/+OHT treatment.

After pre-processing, alignment and sorting, it uses converts the paired end alignment to a single read that starts at the first base of the upstream read and ends at the last base of the downstream read. 

This is then used to calculate read coverage around the cut asisi recognition sites.

This read coverage is normalised to the total read number to give reads per million for each nucleotide specified by the bed file.

Log2 fold change of +OHT(cut)/-OHT(uncut) can also be calculated to show the "damage-induced coverage".

The read coverage files can then be passed to the diva_plotting_pipe.R Rscript that generates a variety of plots.

# Plotting with diva_plotting_pipe.R

Dependencies:
dplyr
zoo
ggplot2
gplots

This script takes three positional arguments:
1. The directory in which the coverage file/bed file are saved.
2. The name of the coverage file.
3. A custom bed file of the sites to plot with the following columns: chromosome, start coodinate of cut site, end coordinate of cut site, sites names(normally just numbered), relative transcriptional activity.

All plots generated will be saved in the directory of the coverage file.

The following plots are generated:
1. Metagenes of a 10kb region centred on the break sites.
2. Boxplots of the same data as the metagenes are plotted by averaging the read coverage in a 1kb region centred on the break sites.
3. A correlation plot of coverage(same data as the boxplots) against relative transcriptional activity.
4. A heatmap of a 40kb region centred on the break sites and ordered from highest transcriptional activity to lowest. 

In addition, csv files are written continaing the data used to generate the boxplots, the associated statistical tests and the correlation tests.
