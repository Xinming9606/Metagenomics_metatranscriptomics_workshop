#!/usr/bin/bash

## QC metatranscriptomic data, split into rRNA, mRNA

## default settings

java_mem=24
threads=8

## refs

refPhix="/home/ubuntu/workshop/refs/phiX174_virus.fa"
refRNA="/home/ubuntu/workshop/refs/silva_rfam_all_rRNAs.fa"
refAdp="/home/ubuntu/workshop/refs/adapters.fa"

# Activate env

eval "$(conda shell.bash hook)"
conda activate /home/ubuntu/conda/bbtools

## mark the start time

starTime=`date +%s`

wkdir=$1
outDir=$2

sample="SRR12632504"

## S3, split into rRNA, mRNA

mkdir -p ${outDir}/03_split
mkdir -p ${outDir}/03_split/mRNA
mkdir -p ${outDir}/03_split/rRNA
mkdir -p ${outDir}/03_split/PhiX

bbsplit.sh \
    in1=${outDir}/02_filt/${sample}_R1.fastq.gz \
    in2=${outDir}/02_filt/${sample}_R2.fastq.gz \
    ref_rRNA=${refRNA} \
    ref_phix=${refPhix} \
    basename="${outDir}/03_split/${sample}_%_R#.fastq.gz" \
    outu1=${outDir}/03_split/mRNA/${sample}_R1.fastq.gz \
    outu2=${outDir}/03_split/mRNA/${sample}_R2.fastq.gz \
    maxindel=20 minratio=0.65 minhits=1 \
    refstats=${outDir}/logs/stats_split_${sample}.txt \
    threads=${threads} \
    k=13 local=t machineout=t \
    pigz=t unpigz=t ziplevel=9 \
    -Xmx${java_mem}G \
    1> ${outDir}/logs/log_03_split_${sample}.log \
    2> ${outDir}/logs/log_03_split_${sample}.err

    mv ${outDir}/03_split/*_phix*.fastq.gz ${outDir}/03_split/PhiX
    mv ${outDir}/03_split/*_rRNA*.fastq.gz ${outDir}/03_split/rRNA

### Report the running time

conda deactivate

## mark the end time

endTime=`date +%s`

## running time summary

runTime=$( echo "($endTime - $starTime + 59)/60" | bc -l )

echo -e "The flow duration: ${runTime%%.*} minutes."
