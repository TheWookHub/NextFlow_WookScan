#!/usr/bin/env python
import phippery as ph
import argparse

parser = argparse.ArgumentParser()

# Making use of the ".phip" file in pickle_data from Phip-flow output
parser.add_argument("-i", help = ".phip file of phip-flow output.", type=str)
parser.add_argument("-u_pep_id", help = "user defined peptide id tag.", type=str)
args = parser.parse_args()

# we load the data here and then take out edgeR hits, also making some pandas dataframe
phip_result = ph.load(args.i)
# phip_sample_table = phip_result.sample_table.to_pandas()
phip_edgeRhits = phip_result.edgeR_hits.to_pandas()
# phip_edgeRhits = phip_result.edgeR_hits.to_pandas().rename(columns = phip_sample_table.technical_replicate_id)

# here we reset the index so it becomes back to one of the columns
# we can use this later to merge data when matching IDs.
phip_edgeRhits.reset_index(inplace = True)

# The beads only samples are give as NaN, so we will put everything to False
# Note - AVARDA works with both 1/0 or True/False. You can change to one that
# you prefer. 
phip_edgeRhits.fillna(False, inplace = True)

# here we extract the peptide annotate information. Main focus it to match original_id / phippery-happy id and the u_pep_id
# which we are going to construct
phip_anno_table = ph.get_annotation_table(phip_result, dim = 'peptide')

# we extract the relevent columns for our id matching:
if('original_id' in phip_anno_table):
    wookscan_table = phip_anno_table[['original_id','oligo']].reset_index()
    # and we get something like this (note peptide_metadata is the index, not treated as a column):
    #
    #       peptide_metadata  peptide_id  original_id                                              oligo
    #       0                          0            1  ATGCGCAGCTTGCTGTTTGTGGTCGGTGCTTGGGTCGCTGCTCTCG...
    #       1                          1            2  ACTACAACCACCGCTGCCGCAGGGAACACATCTGCAACAGCTTCTC...
    #       2                          2            3  ATTACCGCTGCCGCTCCTCCAGGTCATTCAACACCTTGGCCTGCAC...    
    #       ...                      ...          ...                                                ...    
    #       128255                128255       128286  ATCCCTGCCAGCAACGAAACGGATAATAGCCCACTGGGGGGGTATA...
    #       128256                128256       128287  ATTATGACAAGTTCAAAATTTGGCGGGGTCAATGTTTGGAATCGCT...
    
    # we create the new u_pep_id columsn in wookscan_table
    # 'WOOKSCAN_001_' is up to your analysis tag. You can 
    # name it 'UNICORN_069_' if it suits your need.
    wookscan_table['u_pep_id'] = [args.u_pep_id +"_" + x for x in wookscan_table.oligo.values]

    # our AVARDA ready table will be made like via merging:
    phip_edgeRhits_ready = wookscan_table.loc[
        :,
        ['peptide_id','original_id','u_pep_id']
    ].merge(
        phip_edgeRhits, 
        on = 'peptide_id'
    ).drop(
        ['peptide_id','original_id'],
        axis = 1
    ).drop_duplicates()

    # Create virlib Table from the phippery output for AVARDA
    virlib_table = wookscan_table.loc[:,['u_pep_id','original_id']].rename(columns = {'original_id':'pep_id'})

else:
    wookscan_table = phip_anno_table[['oligo']].reset_index()

    # or if there's no 'original_id' in the columns then:
    #
    #       peptide_metadata  peptide_id                                              oligo
    #       0                          0  aggaattctacgctgagtATGTTCCTGATCCTGCTGATCTCTCTGC...
    #       1                          1  aggaattctacgctgagtTGCACCCTGGACCCGCGTCTGAAAGGTT...
    #       2                          2  aggaattctacgctgagtATCTCTACCGACACCGTTGACGTTACCA...    
    #       ...                      ...                                                ...    
    #       10045                  10045  aggaattctacgctgagtACCAACGACCCGATCCGTTTCTGCCTGG...
    #       10046                  10046  aggaattctacgctgagtACCGTTTGCAAAGTTTGCGGTTGCTGGC...

    # we create the new u_pep_id columsn in wookscan_table
    # 'WOOKSCAN_001_' is up to your analysis tag. You can 
    # name it 'UNICORN_069_' if it suits your need.
    wookscan_table['u_pep_id'] = [args.u_pep_id +"_" + x for x in wookscan_table.oligo.values]

    # our AVARDA ready table will be made like via merging:
    phip_edgeRhits_ready = wookscan_table.loc[
        :,
        ['peptide_id','u_pep_id']
    ].merge(
        phip_edgeRhits, 
        on = 'peptide_id'
    ).drop(
        ['peptide_id'],
        axis = 1
    ).drop_duplicates()

    # Create virlib Table from the phippery output for AVARDA
    virlib_table = wookscan_table.loc[:,['u_pep_id','peptide_id']].rename(columns = {'peptide_id':'pep_id'})


# We now write both files to file
phip_edgeRhits_ready.to_csv(
    "PhipperyEdgeRHITS_AVARDA_Input.csv",
    header = True,
    index = False
)

virlib_table[
    ~virlib_table.u_pep_id.duplicated()
].to_csv(
    args.u_pep_id + "_virlib_names.csv",
    index = False,
    header = True
)
