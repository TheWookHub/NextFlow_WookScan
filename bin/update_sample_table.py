#!/usr/bin/env python

import argparse,re,os
import pandas as pd

# parsed arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", type=str) # original sample table
parser.add_argument("-t", type=str) # filtered sample names
parser.add_argument("-r", type=str) # results directory of where filtered fastq are copied
parser.add_argument("-o", type=str) # output name
args = parser.parse_args()

# nextflow config params for where the filtered fastq are located
FILTERED_RESULT_PATH = args.r

# read in original sample table
og_sample = pd.read_csv(args.s, header=0)
if "fastq_filepath" not in og_sample.columns:
    raise KeyError("Must provide 'fastq_filepath' column in the table")

# read in filtered sample names
filtered = pd.read_csv(    
    args.t,
    header = None
).rename(
    columns = {0:'Location'}
)

# extract the file base names from filtered and original sample tables
# and calling them as "Base". We're also changing the file paths to point to filtered
# fastq files in the results directory.
filtered['Base'] = [re.split('_filtered',re.split('/',x)[-1])[0] for x in filtered.Location]
filtered['Results_Location'] = [
    os.path.join(FILTERED_RESULT_PATH,x+"_filtered.fastq.gz") for x in filtered.Base
]
og_sample['Base'] = [re.split('\\.',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath]

# collect the other columns from original sample table that's not
# "Base" or "fastq_filepath"
other_cols = [x for x in og_sample.columns.to_list() if x not in ['Base','fastq_filepath']]

# merge on "Base" column
merged_df = pd.merge(og_sample, filtered, on='Base')

# Get the new file location and other columns then
# rename the column to "fastq_filepath" so that we can output
# a new filtered sample table pointing to the correct filtered fastqs
filtered_final = merged_df.loc[:,['Results_Location'] + other_cols]
filtered_final.rename(columns={'Results_Location':'fastq_filepath'},inplace=True)
filtered_final.to_csv(args.o, index=False, na_rep="NA")