#!/usr/bin/python3

import os
import pandas as pd
import glob
import argparse


## functions

def read_pileup_coverage(pileup_file, coverage_measure="Median_fold"):
    filename = os.path.basename(pileup_file)
    stem, ext = os.path.splitext(filename)
    sampleName = stem.replace('coverage_', '')
    df = pd.read_csv(pileup_file, index_col=0, sep="\t",
         usecols=["#ID", coverage_measure, "Plus_reads", "Minus_reads"]
         )
    df.index.names = ["Contig"]
    df.loc[df[coverage_measure] < 0, coverage_measure] = 0
    df.insert(1, 'genome', df.index.str.split(pat='_',n=1).str[0])
    df.eval("Reads = Plus_reads + Minus_reads", inplace=True)
    df.drop(["Plus_reads", "Minus_reads","Median_fold"], axis=1, inplace=True)
    df2 = df.groupby('genome').sum()
    df2.insert(0, 'Sample', sampleName)
    return df2

def main():
    parser = argparse.ArgumentParser(
        description="get the coverage of mRNA reads mapped to genomes over different samples"
    )

    parser.add_argument("-f", "--filePathPattern", required=True,
            help="file path pattern where the pileup coverage table generated from BAM/SAM files")
    parser.add_argument("-o", "--output", required=True, help="The output file name")

    args = parser.parse_args()

    # Load data
    
    files = glob.glob(args.filePathPattern)

    files.sort() ## sort

    ## read in table and combine

    dfList = list()

    for f in files:
        df = read_pileup_coverage(f, coverage_measure="Median_fold")
        dfList.append(df)

    df = pd.concat(dfList, axis=0).reset_index()
    df_wide = df.pivot(index='genome', columns='Sample', values='Reads')
    print(df_wide.shape)
    print(df_wide.head())

    df_wide.to_csv(args.output, sep='\t', index = True)

if __name__ == "__main__":
    main()





