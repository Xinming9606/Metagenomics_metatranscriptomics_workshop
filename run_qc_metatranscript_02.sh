#!/usr/bin/bash

## QC metatranscriptomic data, filter adapters, low quality, short length ones

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

## S2, filter adapters, low quality, short length ones
    
mkdir -p ${outDir}/02_filt
   
bbduk.sh \
    in1=${outDir}/01_refmt/${sample}_R1.fastq.gz \
    in2=${outDir}/01_refmt/${sample}_R2.fastq.gz \
    out1=${outDir}/02_filt/${sample}_R1.fastq.gz \
    out2=${outDir}/02_filt/${sample}_R2.fastq.gz \
    ref=${refAdp} \
    interleaved=f \
    stats=${outDir}/logs/stats_filt_${sample}.txt \
    overwrite=true \
    trd=t hdist=1 k=23 ktrim=r \
    mink=11 trimq=10 qtrim="r" \
    threads=${threads} \
    minlength=50 \
    maxns=0 minbasefrequency=0 \
    ecco=t prealloc=t \
    pigz=t unpigz=t \
    -Xmx${java_mem}G \
    2> ${outDir}/logs/log_02_filt_${sample}.err \
    1> ${outDir}/logs/log_02_filt_${sample}.log

### Report the running time

conda deactivate

## mark the end time

endTime=`date +%s`

## running time summary

runTime=$( echo "($endTime - $starTime + 59)/60" | bc -l )

echo -e "The flow duration: ${runTime%%.*} minutes."
