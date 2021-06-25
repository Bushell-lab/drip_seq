2#!/bin/bash
#2>&1 | tee

srcdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" # srcdir is the location of this script and should contain all other scripts
rawdir=$PWD # rawdir is the location the script was run from and should only contain the bcl directory
samples=("$@") # arguments are the list of samples (not marked for damaged or undamaged)

undmg=()
dmg=()
for ((i=0;i<${#samples[@]};i++))
do
    undmged=$(echo ${samples[i]}0)
    undmg[i]=$undmged
    dmged=$(echo ${samples[i]}4)
    dmg[i]=$dmged
done
wait
sample_ID=("${undmg[@]}" "${dmg[@]}")

################################################################# Initial setup and fastq unpacking ####################################################################

echo -e Setting up directory structure"\n"

mkdir $rawdir/log
mkdir $rawdir/processed
mkdir $rawdir/rawdata
mkdir $rawdir/rawdata/fastq
mkdir $rawdir/rawdata/fastqc
mkdir $rawdir/rawdata/pe_int
mkdir $rawdir/processed/covs
mkdir $rawdir/processed/wigs





echo -e "\n"Input bcl directory:
read bcldir # user inputs the bcl folder
echo -e Generating fastqs"\n"
bcl2fastq -R $bcldir -o $rawdir/rawdata/fastq -r 15 -w 15 -p 30 --no-lane-splitting --ignore-missing-bcls
rm $rawdir/rawdata/fastq/Undetermined_* # remove the undetermined fastq files
echo -e Unzipping fastqs"\n"
for file in $rawdir/rawdata/fastq/*.gz # simultaneously unzip all fastq fiiles
do
  gzip -d ${file} & done
wait










################################################################ QC and trimming ####################################################################

echo -e Removing low quality reads and QC checking"\n"
for sample in ${sample_ID[@]}
do
   fastp -g -w 5 -l 5 -Q -A -h $rawdir/rawdata/fastqc/${sample}_fastp.html -i $rawdir/rawdata/fastq/${sample}_S*_R1_001.fastq -I $rawdir/rawdata/fastq/${sample}_S*_R2_001.fastq --out1 $rawdir/rawdata/fastq/${sample}_qc1.fq --out2 $rawdir/rawdata/fastq/${sample}_qc2.fq & done
wait





########################################################## Alignment and sorting ####################################################################

echo -e Alignment results"\n\n" > $rawdir/log/alignments.log
echo -e "\n\n\n"Alignment results"\n"
for sample in ${sample_ID[@]} # align fastqs and pipe straight to samtools sort, outputting the bowtie alignment results to a log file
do
   (bowtie2 -t -p 60 --fr --maxins 5000 --dovetail -x /mnt/data/R11/raw_sequencing_data/Aldo_raw_seq/genomes/hg38/hg38 \
   -1 $rawdir/rawdata/fastq/${sample}_qc21.fq -2 $rawdir/rawdata/fastq/${sample}_qc22.fq \
   | samtools sort -o $rawdir/rawdata/${sample}sort.bam -O bam -l 0 -@ 60 -m 1G) \
   2> $rawdir/log/${sample}_align.log

   echo -e ${sample} aligned successfully
   echo -e "\n"${sample} alignment: >> $rawdir/log/alignments.log
   tail -n 19 $rawdir/log/${sample}_align.log >> $rawdir/log/alignments.log
done
wait
echo -e "\n"Indexing
for sample in ${sample_ID[@]} # simultaneously index all alignments
do
   samtools index $rawdir/rawdata/${sample}sort.bam & done
wait










echo -e "\n""\n"Sorting by name"\n"
for sample in ${sample_ID[@]}
do
   samtools sort $rawdir/rawdata/${sample}sort.bam -o $rawdir/rawdata/pe_int/${sample}nsort.bam -n -O bam -l 0 -@ 60 -m 1G
   echo -e ${sample} sorted
done
wait
echo -e "\n""\n"Fixing samples"\n"
for sample in ${sample_ID[@]}
do
   samtools fixmate $rawdir/rawdata/pe_int/${sample}nsort.bam $rawdir/rawdata/pe_int/${sample}fixed.bam -@ 4 & done
wait
echo -e converting to bedpe"\n"
for sample in ${sample_ID[@]}
do
   bedtools bamtobed -bedpe -i $rawdir/rawdata/pe_int/${sample}fixed.bam > $rawdir/rawdata/pe_int/${sample}.bedpe & done
wait
echo -e cutting bedpe to bed"\n"
for sample in ${sample_ID[@]}
do
    bedpeTobed.py $rawdir/rawdata/pe_int/${sample}.bedpe $rawdir/rawdata/pe_int/${sample}cut.bed & done
wait
echo -e converting to bam"\n"
for sample in ${sample_ID[@]}
do
   bedtools bedtobam -i $rawdir/rawdata/pe_int/${sample}cut.bed -g /mnt/data/R11/raw_sequencing_data/Aldo_raw_seq/genomes/hg38/GRCh38.p10.genome > $rawdir/rawdata/pe_int/${sample}cut.bam & done
wait
echo -e "\n""\n"Sorting samples"\n"
for sample in ${sample_ID[@]}
do
   samtools sort $rawdir/rawdata/pe_int/${sample}cut.bam -l 0 -o $rawdir/rawdata/${sample}_cutsort.bam -O bam -@ 60 -m 1G
done
wait
echo -e "\n""\n"Indexing"\n"
for sample in ${sample_ID[@]}
do
    samtools index $rawdir/rawdata/${sample}_cutsort.bam & done
wait






echo -e calculating coverage"\n"
for sample in ${sample_ID[@]}
do
    samtools depth -aa -d 0 -b $srcdir/5kbreaksites.bed $rawdir/rawdata/${sample}_cutsort.bam > $rawdir/rawdata/${sample}_pe.cov & done
wait
echo -e calculating readcounts"\n"
>$rawdir/log/readcounts.log
for sample in ${sample_ID[@]}
do
   readcount=$(samtools view -c $rawdir/rawdata/${sample}_cutsort.bam)
   echo -e ${sample} has $readcount reads
   echo -e $readcount>>$rawdir/log/readcounts.log
done
readarray -t readcounts < $rawdir/log/readcounts.log

echo -e "\n""\n""\n"normalising covs"\n"
for ((i=0;i<${#sample_ID[@]};i++))
do
    awk -v x=${readcounts[$i]} '{print $1, $2, (($3+1)/(x/1000000))}' $rawdir/rawdata/${sample_ID[$i]}_pe.cov > $rawdir/rawdata/${sample_ID[$i]}_penorm.cov & done
wait

echo -e logging samples"\n"
for ((i=0;i<${#dmg[@]};i++))
do
    awk 'NR==FNR{a[FNR]=$3; next} {print $1, "\t", $2, "\t", log($3/a[FNR])/log(2)}' $rawdir/rawdata/${undmg[$i]}_penorm.cov $rawdir/rawdata/${dmg[$i]}_penorm.cov > $rawdir/processed/covs/${samples[$i]}_pe_log.cov & done
wait

