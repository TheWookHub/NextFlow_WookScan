#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser()

# Making use of the ".phip" file in pickle_data from Phip-flow output
parser.add_argument("-i", help = "sample table used for phippery.", type=str)
parser.add_argument("-p", help = "output path.", type=str)
args = parser.parse_args()

# ---------- Function(s) ---------- #
# Linking sample files (fastq) back to the subject names of
# AVARDA output
def avardaSampleID_Table(sample_table_file,output_prefix):
    sample_table = pd.read_table(sample_table_file, sep = ',')
    st = sample_table.loc[
        :,['fastq_filepath']
    ].reset_index().rename(
        columns = {'index':'name'}
    )
    st['name'] = ['X' + str(x) for x in st.name]
    fileName = output_prefix + "_runSampleID_Table.csv"
    st.to_csv(fileName, index = False, header = True)

# ---------- Main ---------- #

avardaSampleID_Table(args.i,args.p)