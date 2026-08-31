#!/usr/bin/env python

import argparse,re,os
import pandas as pd

# parsed arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", type=str) # original sample table
parser.add_argument("-t", type=str) # filtered sample names
parser.add_argument("-r", type=str) # results directory of where filtered fastq are copied
parser.add_argument("-o", type=str) # output name
parser.add_argument("-m", type=str, default="trimmed", choices=["trimmed","filtered","assembled"]) # declare if we're at trimming stage or filtering stage.
args = parser.parse_args()

# nextflow config params for where the filtered fastq are located
RESULT_PATH = args.r

# identifying naming for trimmed or filtered
stage = args.m

# read in original sample table
og_sample = pd.read_csv(args.s, header=0)
if(stage in ["trimmed","filtered"]):
    if ("fastq_filepath" not in og_sample.columns and "fastq_filepath2" not in og_sample.columns):
        raise KeyError("Must provide 'fastq_filepath' and 'fastq_filepath2' columns in the table")
    # read in sample names
    trimmed_or_filtered = pd.read_csv(    
        args.t,
        header = None,
        sep = "\t"
    ).rename(
        columns = {0:'R1_Path', 1:'R2_Path'}
    )
    print(trimmed_or_filtered)
    # extract the file base names from filtered and original sample tables
    # and calling them as "Base". We're also changing the file paths to point to filtered
    # fastq files in the results directory.
    trimmed_or_filtered['R1_Base'] = [re.split(f'_{stage}',re.split('/',x)[-1])[0] for x in trimmed_or_filtered.R1_Path]
    trimmed_or_filtered['R2_Base'] = [re.split(f'_{stage}',re.split('/',x)[-1])[0] for x in trimmed_or_filtered.R2_Path]

    trimmed_or_filtered['R1_Results_Location'] = [
        os.path.join(RESULT_PATH,x+f"_{stage}.fastq.gz") for x in trimmed_or_filtered.R1_Base
    ]
    trimmed_or_filtered['R2_Results_Location'] = [
        os.path.join(RESULT_PATH,x+f"_{stage}.fastq.gz") for x in trimmed_or_filtered.R2_Base
    ]

    if(stage == "trimmed"):
        og_sample['R1_Base'] = [re.split('\\.',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath]
        og_sample['R2_Base'] = [re.split('\\.',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath]
    elif(stage == "filtered"):
        og_sample['R1_Base'] = [re.split('_trimmed',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath]
        og_sample['R2_Base'] = [re.split('_trimmed',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath2]

    # collect the other columns from original sample table that's not
    # "Base" or "fastq_filepath"
    other_cols = [x for x in og_sample.columns.to_list() if x not in ['R1_Base','fastq_filepath','R2_Base','fastq_filepath2']]

    # merge on "Base" column
    merged_df = pd.merge(og_sample, trimmed_or_filtered, left_on='R1_Base', right_on='R1_Base')

    # Get the new file location and other columns then
    # rename the column to "fastq_filepath" so that we can output
    # a new filtered sample table pointing to the correct filtered fastqs
    final = merged_df.loc[:,['R1_Results_Location','R2_Results_Location'] + other_cols]
    final.rename(columns={'R1_Results_Location':'fastq_filepath', 'R2_Results_Location':'fastq_filepath2'}, inplace=True)    

elif(stage == "assembled"): # reason we need this is because at this stage, we've assembled PE->SE. 
    if ("fastq_filepath" not in og_sample.columns):        
        raise KeyError("Must provide 'fastq_filepath' column in the table")    
    
    og_sample.drop(columns=['fastq_filepath2'], inplace=True)
    
    assembled = pd.read_csv(    
        args.t,
        header = None
    ).rename(
        columns = {0:'Location'}
    )
    assembled['Base'] = [re.split(f'.{stage}',re.split('/',x)[-1])[0] for x in assembled.Location]
    assembled['Results_Location'] = [
        os.path.join(RESULT_PATH,x+f".{stage}.fastq.gz") for x in assembled.Base
    ]
    og_sample['Base'] = [re.split('_R1_.',re.split('/',x)[-1])[0] for x in og_sample.fastq_filepath]
    
    # collect the other columns from original sample table that's not
    # "Base" or "fastq_filepath"
    other_cols = [x for x in og_sample.columns.to_list() if x not in ['Base','fastq_filepath']]

    # merge on "Base" column
    merged_df = pd.merge(og_sample, assembled, on='Base')

    # Get the new file location and other columns then
    # rename the column to "fastq_filepath" so that we can output
    # a new filtered sample table pointing to the correct filtered fastqs
    final = merged_df.loc[:,['Results_Location'] + other_cols]
    final.rename(columns={'Results_Location':'fastq_filepath'},inplace=True)
final.to_csv(args.o, index=False, na_rep="NA")

