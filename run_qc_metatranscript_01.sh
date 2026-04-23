#!/usr/bin/bash

## QC metatranscriptomic data, reformat data

## default settings

java_mem=24
threads=8


refPhix="/home/ubuntu/workshop/refs/phiX174_virus.fa"
refRNA="/home/ubuntu/workshop/refs/silva_rfam_all_rRNAs.fa"
refAdp="/home/ubuntu/workshop/refs/adapters.fa"

eval "$(conda shell.bash hook)"

conda activate /home/ubuntu/conda/bbtools

## mark the start time

starTime=`date +%s`

wkdir=$1
outDir=$2

sample="SRR12632504"

fp_raw="/home/ubuntu/workshop/metatranscript/raw"

# S1, reformat input

mkdir -p ${outDir}/logs
mkdir -p ${outDir}/01_refmt

reformat.sh \
    in1=${fp_raw}/${sample}_R1.fastq.gz \
    in2=${fp_raw}/${sample}_R2.fastq.gz \
    out1=${outDir}/01_refmt/${sample}_R1.fastq.gz \
    out2=${outDir}/01_refmt/${sample}_R2.fastq.gz \
    interleaved=f \
    qtrim=r trimq=10 minlength=50 \
    overwrite=true \
    verifypaired=t \
    threads=${threads} \
    -Xmx${java_mem}G \
    2> ${outDir}/logs/log_01_rfmt_${sample}.txt

conda deactivate

### Report the running time

## mark the end time

endTime=`date +%s`

## running time summary

runTime=$( echo "($endTime - $starTime + 59)/60" | bc -l )

echo -e "The flow duration: ${runTime%%.*} minutes."
