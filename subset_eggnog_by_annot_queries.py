#!/usr/bin/python3

## This script subset eggnog anntation table based on queried list of KOs
## The counts can be linked and even convert to TPM or RPKM
## Sizhong Yang, ZALF, 11.04.2026


import re, argparse, os
import pandas as pd


### some functions

def normalize_counts(df, gene_len, method="TPM"):
    # Work on a copy to avoid modifying original data
    df = df.copy()
    gene_len = gene_len.copy()
    # Ensure GeneID is a column
    if df.index.name == 'GeneID':
        df = df.reset_index()
    # Standardize gene_length column
    gene_len = gene_len.rename(columns={ 'geneLen': 'gene_length', 'gene_len': 'gene_length'})
    # Merge
    newDF = pd.merge(df, gene_len, on='GeneID', how='left')
    # Identify numeric columns
    numCol = [x for x in newDF.columns if x not in ['GeneID', 'gene_length']]
    # Step 1: normalize by gene length (kb)
    A = newDF.loc[:, numCol].div(newDF.gene_length, axis=0).multiply(1e3)
    if method == "TPM":
        # TPM: normalize after length scaling
        scaling = A.sum(axis=0, numeric_only=True)
        result = A.div(scaling, axis=1).multiply(1e6)
    elif method == "RPKM":
        # RPKM: normalize by library size
        total_reads = newDF.loc[:, numCol].sum(axis=0, numeric_only=True)
        result = A.div(total_reads, axis=1).multiply(1e6)
    else:
        raise ValueError("method must be 'TPM' or 'RPKM'")
    # Add GeneID back
    result['GeneID'] = newDF['GeneID']
    return result

def flatten(df, query_col, query_list, sep=","):
    # Convert query list to set (O(1) lookup)
    query_set = set(query_list)
    # Split into lists
    df = df.copy()
    df[query_col] = df[query_col].str.split(sep)
    # Explode into rows
    df = df.explode(query_col)
    # Filter matching values
    df = df[df[query_col].isin(query_set)]
    return df

## Work on taxonomy and collapse data

parser = argparse.ArgumentParser(description='Extract the eggNOG annoation result by queried items.')
parser.add_argument('--eggnog','-e', help='the input table of eggNOG annotation', required = True)
parser.add_argument("--queries","-q",help='the list of queries items', required = True) 
parser.add_argument('--output','-o', help='[optional] write the table out to file', type=str, required = False)
parser.add_argument('--collaps_output','-m', help='[optional] write to the file where the collapsed table by queried items', type=str, required = False)
parser.add_argument("--counts", "-c", help="the input table of the raw counts", type = str, required = False)
parser.add_argument("--genelen",'-g', help="the gene_length table, which must be provided if normalize to TPM or RPKM", type = str, required = False)
parser.add_argument("--norm", "-n", choices=["TPM", "RPKM", "none"], default="none", type= str,
     help="Normalization method: TPM, RPKM, or none (default: none)")
parser.add_argument("--sep", "-s", default="\t", type= str,
     help="seperator of the data frame (default: TAB '\t')")

### parse arguments

args = parser.parse_args()

### Processing from here
## read in eggnog table

select_cols = ['Query', 'COG_cat', 'Description', 'Name',
        'EC', 'KO', 'KEGG_Module', 'PFAMs']
new_colname = {'Query':'GeneID',
        'COG_cat':'COG_category',
        'KOs':'KO', 'KEGG_KOs':'KO',
        'Description':'function_descrip',
        'Name': 'gene',}

eggnog = pd.read_csv(args.eggnog, header =0, sep = args.sep, usecols = select_cols, comment='#', dtype=str)

eggnog.rename(columns = new_colname, inplace =True)

print('\nThe dimension of the input eggnog file:', eggnog.shape)

eggnog = eggnog[eggnog.KO.ne('-')]
eggnog['KO'] = eggnog['KO'].str.replace("ko:","")
geneID_list = eggnog['GeneID']

nr,nc = eggnog.shape
print('\nThe eggNog data include %s rows, %s columns with clear annotated KOs'%(nr,nc))

## read in the queries

query_list = list()

with open(args.queries) as f:
    query_list = list({line.strip() for line in f if line.strip()})

print('\nA total of %s KO queries'%(len(query_list)))

### filter by queries querys

df_subset = flatten(df = eggnog, query_col= "KO", query_list=query_list, sep= args.sep)

### including the counts

if args.counts:
    counts = pd.read_csv(args.counts, header=0, sep='\t')
    if counts.columns[0].startswith("Unnamed"): # if first column has no name → pandas calls it "Unnamed: 0"
        counts = counts.rename(columns={counts.columns[0]: "GeneID"})
    else:
        counts.rename(columns = {'GeneNr':'GeneID','geneID':'GeneID','gene_id':'GeneID'}, inplace = True)
    tmp = counts.set_index('GeneID')
    tmp2 = tmp.reindex(geneID_list)
    # normalization 
    if args.norm in ["TPM","RPKM"]:
        print(f"\nNormalization method: {args.norm}\n")
        if not args.genelen:
            print('\nGene_length is essential to calculate TPM, will skip calculate TPM but retains raw counts')
            subCount = tmp2
        else:
            genel_df = pd.read_csv(args.genelen, header=None, names=['GeneID','gene_length'],sep='\t')
            genel_df.rename(columns = {'GeneNr':'GeneID','geneID':'GeneID','gene_id':'GeneID',
                'geneLen': 'gene_length','gene_len':'gene_length'}, inplace = True)
            subCount = normalize_counts(tmp2, genel_df, method=args.norm)
    else: 
        print("\nNormaliation is not requested. Keep raw counts.\n")
        subCount = tmp2.reset_index()
    if not subCount.empty:
        final_cols = ["KO"] + [col for col in subCount.columns if col != "KO" and col != 'GeneID']
        res = df_subset.merge(subCount, on="GeneID", how="left")  
    else:
        res = df_subset
        final_cols = None
    ## collapse data by the query_col if final_cols
    if final_cols:
        result = res[final_cols].groupby("KO").sum()
        print(result.shape)
        if args.collaps_output:
            print(f"\nFile was written to {args.collaps_output}\n")
            result.to_csv(args.collaps_output, index = True,sep='\t')
    if args.output:
        print(f"File was written to {args.output}\n")
        res.to_csv(args.output, index = False,sep='\t')
    else:
        pass
        #print(res.shape)
        #print(res.head())
else:
    #print(df_subset.shape)
    #print(df_subset.head())
    if args.output:
        print(f"\nFile was written to {args.output}\n")
        df_subset.to_csv(args.output, index = False, sep='\t')
