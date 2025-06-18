#!/usr/bin/env python

import argparse,re
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-s", type=str) # original sample table
parser.add_argument("-t", type=str) # trimmed sample names
parser.add_argument("-o", type=str) # output name
args = parser.parse_args()

# read in original sample table
og_sample = pd.read_csv(args.s, header=0)
if "fastq_filepath" not in og_sample.columns:
    raise KeyError("Must provide 'fastq_filepath' column in the table")

# read in trimmed sample names
trimmed = pd.read_csv(    
    args.t,
    header = None
).rename(
    columns = {0:'Location'}
)

# extract the file base names from trimmed and original sample tables
# and calling them as "Base"
# trimmed['File'] = [re.split('/',x)[-1] for x in trimmed.Location]
# trimmed['Base'] = [re.split('_trimmed',x)[0] for x in trimmed.File]
trimmed['Base'] = [re.split('_trimmed',re.split('/',x)[-1])[0] for x in trimmed.Location]
og_sample['Base'] = [re.split('\\.',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath]

# collect the other columns from original sample table that's not
# "Base" or "fastq_filepath"
other_cols = [x for x in og_sample.columns.to_list() if x not in ['Base','fastq_filepath']]

# merge on "Base" column
merged_df = pd.merge(og_sample, trimmed, on='Base')

# Get the new file location and other columns then
# rename the column to "fastq_filepath" so that we can output
# a new trimmed sample table pointing to the correct trimmed fastqs
trimmed_final = merged_df.loc[:,['Location'] + other_cols]
trimmed_final.rename(columns={'Location':'fastq_filepath'},inplace=True)
trimmed_final.to_csv(args.o, index=False, na_rep="NA")