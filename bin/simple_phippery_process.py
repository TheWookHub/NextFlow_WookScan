#!/usr/bin/env python

import phippery as ph; import pandas as pd; import datatable as dt
import argparse

#####################
# default variables #
#####################

status_control = 'beads_only'
status_observation = 'empirical'
library = 'library'

relevant_cols = [
    'technical_replicate_id','control_status',
    'raw_total_sequences','reads_mapped',
    'percent_mapped','percent_peptides_detected',
    'percent_peptides_between_10_and_100'
]

###########################
# defining some functions #
###########################


def getPeptideSummary(sample_list,hits_counts_df):    
    summaryList = []
    for sample in sample_list:
        r = hits_counts_df.loc[
            :,['Species',sample]
        ].groupby(
            'Species'
        ).sum().sort_values(
            by = sample,
            ascending = False
        )
        summaryList.append(r)
    return pd.concat(summaryList,axis=1)


# Lazy function so that we could easily rename columns of
# the dataframe to be a bit more meaningful than just numbers
def makeColumnDict(species_list, column_list):
    assert(len(species_list) == len(column_list))
    myDict = {}
    for index in range(0,len(species_list)):
        myDict[column_list[index]] = species_list[index]
    return myDict


###########################
# main script starts here #
###########################

# dealing with arguments
parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str) # data.phip file
args = parser.parse_args()

# loading the phip data
phip_data = ph.load(args.d)
peptide_table = phip_data.peptide_table.to_pandas()
sample_table = phip_data.sample_table.to_pandas()

# extracting relevant columns from sample table
sample_table_2 = sample_table.loc[:,relevant_cols]

# getting all the sample names, and filtering out those that are
# background controls
all_samples = sample_table_2.technical_replicate_id.to_list()

observation_only = sample_table_2[
    sample_table_2.control_status == status_observation
].technical_replicate_id.to_list()

control_only = sample_table_2[
    sample_table_2.control_status == status_control
].technical_replicate_id.to_list()

reNameDict = makeColumnDict(all_samples,phip_data.counts.to_pandas().columns.tolist())

# retrieve hits table (a 1 or 0 table)
hits = phip_data.edgeR_hits.to_pandas().rename(columns = reNameDict).loc[:,observation_only].astype(int)

#annotated hits
anno_hits = peptide_table.loc[
    :,
    ['original_id','Species']
].merge(
    hits,
    left_index=True,
    right_index=True
)
# retrieve annotated counts table (read count table)
anno_counts = peptide_table.loc[
    :,
    ['original_id','Species']
].merge(
    phip_data.counts.to_pandas(),
    left_index=True,
    right_index=True
).rename(columns = reNameDict)

# annoated hits counts 
# counts that were not a hit will be flattened to 0
anno_hits_counts = peptide_table.loc[
    :,['original_id','Species']
].merge(
    anno_counts.drop(
        columns = control_only + ['original_id','Species']
    ).mul(hits), 
    left_index = True,
    right_index = True
)

# Summarise how many hits counts per species
hits_counts_species = getPeptideSummary(observation_only, anno_hits_counts)
# summary how many peptide hits per species
peptide_hit_per_species = anno_hits.loc[
    :,['Species'] + observation_only
].groupby('Species').sum()

# write the outputs 
sample_table_2.to_csv("sample_table.csv")
anno_hits.to_csv("hits.csv")
anno_hits_counts.to_csv("hits_counts.csv")
hits_counts_species.to_csv("hits_counts_species.csv")
peptide_hit_per_species.to_csv("peptide_hit_per_species.csv")