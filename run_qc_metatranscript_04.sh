#!/usr/bin/bash

## QC metatranscriptomic data, dereplicate mRNA to get uniq seqs for assembly

## default settings

java_mem=24
threads=8

## Activate env

eval "$(conda shell.bash hook)"
conda activate /home/ubuntu/conda/bbtools

## mark the start time

starTime=`date +%s`

wkdir=$1
outDir=$2

sample="SRR12632504"
 
## S4,dereplicate, uniq seqs for assembly, clean seqs for coverage

mkdir -p ${outDir}/04_derep

clumpify.sh \
    in1=${outDir}/03_split/mRNA/${sample}_R1.fastq.gz \
    in2=${outDir}/03_split/mRNA/${sample}_R2.fastq.gz \
    out1=${outDir}/04_derep/${sample}_R1.fastq.gz \
    out2=${outDir}/04_derep/${sample}_R2.fastq.gz \
    overwrite=true \
    dedupe=t \
    threads=${threads} \
    pigz=t unpigz=t \
    -Xmx${java_mem}G \
    2> ${outDir}/logs/log_04_drep_${sample}.err;

### Report the running time

conda deactivate

## mark the end time

endTime=`date +%s`

## running time summary

runTime=$( echo "($endTime - $starTime + 59)/60" | bc -l )

echo -e "\nThe flow duration: ${runTime%%.*} minutes.\n"
