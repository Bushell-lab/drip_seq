#!/bin/bash

# enter the directory which contains the following directory tree and files:
# The directory should have the following subdirectories: rawdata, processed, log and source.
# This script, along with the other scripts and required bed files should be in source.
# rawdata should have a subdirectory named fastq which contains the fastq files. 
# logs of the alignments are saved in the log directory.
# final coverage data files and all generated plots are saved in the processed direcotry. 
# The directory tree enables easy navigation of the data and all the intermediate files generated. 

# Here I have used the sample names drip0 and drip4 as the fastq files were named drip0.fastq and drip4.fastq.
# These should be paired, e.g. if you had another set of samples named chip0 and chip4 then undmg=(drip0 chip0) and dmg=(drip4 chip4)
rawdir=
un_dmg=(drip0)
dmg=(drip4)
sample_ID=("${un_dmg[@]}" "${dmg[@]}")

# The parameters here are relatively uninteresting.
# A qscore cutoff of 20 is pretty standard and the length cutoff of 5 is purely to removed rare 1-2nt sequences output by NGS systems that cause problems in the alignemnt.
# The original fastq is retained as this retains the original data. Later on, BAM/SAM files can be deleted to save storage space as they can be re-generated. 
echo -e Trimming unwanted reads"\n"
for sample in ${sample_ID[@]}
do
   cutadapt -q 20 -m 5 -o $rawdir/rawdata/fastq/${sample}qc.fq $rawdir/rawdata/fastq/${sample}.fastq & done
wait


# WARNING I have used the argument -p 31 here as I assigned 31 threads to the alignemnts. 
# WARNING You should change this to the number of threads your computer can spare, one alignment is done at a time. 
# Aligns each sample printing to the terminal when each is complete.
# In addition, the full log of each alignment to written to a log, and the summary for each alignment is written to the alignments.log file.
echo -e Alignment results"\n\n" > $rawdir/log/alignments.log
echo -e "\n\n\n"Alignment results"\n"
for sample in ${sample_ID[@]}
do
   (bowtie2 -p 31 -x /mnt/data/R11/raw_sequencing_data/Aldo_raw_seq/genomes/hg38/hg38 -U $rawdir/rawdata/${sample}qc.fq -S $rawdir/rawdata/${sample}.sam) 2> $rawdir/log/${sample}_align.log
   echo ${sample} aligned
   echo -e "\n"${sample} alignment: >> $rawdir/log/alignments.log
   tail -n 5 $rawdir/log/${sample}_align.log >> $rawdir/log/alignments.log
done
wait 


# This uses samtools to sort and index the aligments. 
# WARNING sort also has thread (-@) and memory (-m) arguments. 
# WARNING the memory is allocated per thread, i.e. if you have 5 threads to spare and 10Gb of RAM, set it to -@ 5 -m 2G
echo -e "\n""\n"Sorting samples"\n"
for sample in ${sample_ID[@]}
do
   samtools sort $rawdir/rawdata/${sample}.sam -l 0 -o $rawdir/rawdata/${sample}sort.bam -O bam -@ 3 -m 1500M & done
wait
echo -e "\n""\n""\n"Indexing"\n"
for sample in ${sample_ID[@]}
do
    samtools index $rawdir/rawdata/${sample}sort.bam & done
wait


# depth takes the sorted SAM and determines the read coverage at every nucleotide in the bed file you pass to it via -b
# A bed file containing all of the coordinates to generate read coverage at must be provided here. Mine is called asisi20k.bed and is saved in the source directory.
# This bed file is also required later on so if you name it differently or save it somewhere else this must also be changed downstream!
# -aa is required! Without this, positions with 0 read coverage at either end of a region are not reported, which will prevent the downstream analysis.
echo -e calculating coverage"\n"
for sample in ${sample_ID[@]}
do
     samtools depth -aa -d 0 -b $rawdir/asisi20k.bed $rawdir/${sample}sort.bam > $rawdir/${sample}20k.cov & done
wait

# The following code is not neccessary to automate in this way, but I find it easier.
# The for loop fills the readcounts array with the number of reads in each sorted SAM, determined using samtools view -c.
# Sometimes this does not work, but you can manually fill the readcounts array from the readcounts.log
echo -e calculating readcounts"\n"
readcounts=()
for sample in ${sample_ID[@]}
do
   readcount=$(samtools view -c $rawdir/${sample}sort.bam)
   echo -e ${sample} has $readcount reads
   readcounts+=$readcount
done
echo -e readcounts are $readcounts"\n" > $rawdir/log/readcounts.log
# Recording the readcounts for each sample is important in case this step needs re-running.

# Here, the first awk script takes the readconut for each sample and divides the coverage of each nucelotide in that sample by the total readcount/1,000,000.
# This converts the raw coverage to reads per million and is neccessary for any successful comparison between two samples. 
echo -e "\n""\n""\n"normalising gcovs"\n"
for ((i=0;i<${#sample_ID[@]};i++))
do
    awk -v x=${readcounts[$i]} '{print $1, $2, (($3+1)/(x/1000000))}' $rawdir/${sample_ID[$i]}20k.cov > $rawdir/${sample_ID[$i]}_20knorm.cov & done
wait

# Here, a second awkscript takes log2 fold change of the undamaged sample over the damaged sample
# This produced a seperate coverage file with log2 FC. 
echo -e logging samples"\n"
for ((i=0;i<${#un_dmg[@]};i++))
do
    awk 'NR==FNR{a[$1]=$3; next} {print $1, "\t", $2, "\t", log($3/a[$1])}' $rawdir/${un_dmg[$i]}_20knorm.cov $rawdir/${dmg[$i]}_20knorm.cov > $rawdir/${dmg[$i]}_20klog.cov & done
wait



# Assuming the R analysis script is saved in source, this loop passes the coverage file to the R script to generate the neccessary plots.
echo -e passing to Rscript"\n""\n""\n"
for sample in ${dmg[@]}
do
    Rscript $rawdir/source/diva_plotting_pipe.R $rawdir/processed ${sample} & done
wait
echo -e "\n\n\n\n\n"script complete"\n"




