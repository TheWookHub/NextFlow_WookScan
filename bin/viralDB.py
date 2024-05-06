#!/usr/bin/python

from Bio import GenBank
from itertools import chain
from multiprocessing import Pool, cpu_count
import pandas as pd
import datatable as dt
import numpy as np
import subprocess as sp
import viralDB_classes as vDBClass
import re,time, shlex

# ========================= cleanGenBankInfo ========================= #
# This function checks the genbank format for strange data entries.
# On some records, they incorrectly use a single ":". This needs to be
# changed to " :: " for Biopython to be able to read it. Here sed is called
# to do the sub.

def cleanGenBankInfo(fileName):
    tempFile = fileName + ".temp"
    sedcmd1 = f"sed 's/Sequencing Technology: /Sequencing Technology :: /g' {fileName} | "
    sedcmd2 = f"sed 's/Assembly method: /Assembly Method :: /g' | "
    sedcmd3 = f"sed 's/Coverage: /Coverage :: /g' > {tempFile}"
    fullcmd = sedcmd1 + sedcmd2 + sedcmd3
    sp.run(fullcmd, shell=True)    
    return(tempFile)

# ========================= getGenBankInfo ========================= #
# Parse .gb files. Extract the information needed and return a dictionary
# that contains simple organism info and protein info (amino acid sequence
# and the nt sequence which translates into the protein).

# Note: We can assume that the first FEATURE in the .gb format for any
# organism entry will always be "source". This allows us to check if it
# is human, not human, or information not available. If it is not human,
# then we go to the next entry.

def getGenBankInfo(fileName,vhdb_human):
    filteredFile = cleanGenBankInfo(fileName)    
    gb_dict = {}    
    refseq_id_array = list(
        chain.from_iterable(
            [x.split(',') for x in vhdb_human.refseq_id.values.tolist()]
        )
    )
    with open(filteredFile, 'r') as handle:
        for record in GenBank.parse(handle):            
            orgClass = vDBClass.OrganismClass(            
                locus = record.locus,
                version = record.version,
                organism = record.organism,
                full_seq = record.sequence,
                accession = record.accession
            )
            for i in range(0,len(record.features)):
                f_key = record.features[i].key
                match f_key:
                    case "source":                        
                        quals = record.features[i].qualifiers
                        for x in range(len(quals)):
                            q_key = quals[x].key
                            match q_key:
                                case "/host=":
                                    host = re.sub('"','',quals[x].value)
                                    orgClass.addHost(host)
                                case "/db_xref=":
                                    taxon_id = re.sub('"','',quals[x].value.split(':')[1])
                                    orgClass.addTaxonID(taxon_id)
                    case "CDS":
                        protClass = vDBClass.ProteinClass(location = record.features[i].location)
                        quals = record.features[i].qualifiers
                        for x in range(len(quals)):
                            q_key = quals[x].key                         
                            match q_key:
                                case "/product=":
                                    protein = re.sub('"','',quals[x].value)
                                    protClass.addProteinName(protein)
                                case "/protein_id=":
                                    protein_id = re.sub('"','',quals[x].value)
                                    protClass.addProteinID(protein_id)
                                case "/translation=":                                    
                                    protein_seq = re.sub('"','',quals[x].value)
                                    protClass.addProteinSeq(protein_seq)
                        orgClass.appendCDS(protClass)
                        del(protClass)                
                if(orgClass.host is None or 'Homo sapiens' not in orgClass.host):                    
                    if(record.accession[0] in refseq_id_array):                        
                        orgClass.addHost('Homo sapiens')
                    else:
                        break            
            if(orgClass.host is not None and 'Homo sapiens' in orgClass.host):                
                if(orgClass.locus not in gb_dict.keys()):
                    gb_dict[orgClass.accession[0]] = orgClass
                else:
                    print(f"Repeated Accession: {orgClass.accession[0]}")
    handle.close()
    sp.run(f'rm {filteredFile}', shell=True)
    return gb_dict

# ref_friends = getGenBankInfo("sequence_RefSeqPlusNeighbours_new.gb", vhdb)

# ========================= getVHDB_Humans  ========================= #
def getVHDB_humans(fileName):
    VHDB = dt.fread(fileName).to_pandas()
    VHDB_humans = VHDB.loc[
        :,
        [
            'refseq id','host name',
            'virus name','virus tax id',
            'evidence'
        ]
    ][VHDB['host name'] == 'Homo sapiens'].rename(
        columns = {
            'refseq id':'refseq_id', 
            'host name':'host_name',
            'virus name':'virus_name',
            'virus tax id':'virus_tax_id'
        }
    ).reset_index(drop = True).copy()
    return VHDB_humans

# ========================= make pep to species table using bitscore or pident ========================= #
# Here we make a peptide id vs species table with values being the maximum value of bitscore or pident.
# make_pep_species_table_v2 uses taxon id in the sseqid column rather than Accession number.
# In this output, row names will be peptide ids and column names will be species names. The "column name"
# of the peptide ids will be arbitarily labled "V1" just because that's how the AVARDA test files
# provided called it.
         
def make_pep_species_table_v2(colname,bitscoreDF):
    # create a m by n dimension dataframe filled with 0.0
    # where column names are species and row names are peptide ids
    d = pd.DataFrame(
        0.0,
        index = bitscoreDF.qseqid.drop_duplicates().reset_index(drop = True).values,
        columns = bitscoreDF.Species.drop_duplicates().reset_index(drop = True).values
    )
    # now we can can fill in the pepid - Species pair where pident is not 0.0
    for i in range(0,len(bitscoreDF)):
        d.loc[
            bitscoreDF.iloc[i].qseqid,
            bitscoreDF.iloc[i].Species
        ] = bitscoreDF.iloc[i][colname]
    d.index.name = 'V1'
    return d

# ========================= virus to virus peptide set comparison ========================= #
# Conduct a pair-wise peptide set comparison (n^2) between viruses. The table needs to store
# the total probability and unique probability of peptides in the VIR peptide pool. The peptides
# binned to each virus has been bitscore filtered (> 80) to be called 'evidence peptide' for the
# existence of that viral infection.
#
#
# If using WOOKIES_virlib_names_no91273.csv as the VIR library, then the size is 115752 peptides.
# More detailed explanation:
#   Unique Probability  ->  When comparing 2 viruses (e.g. A and B) and their peptide sets, 
#                           there can be three subsets that results. First set is the
#                           intersection between A and B. Second is the difference of A and B.
#                           Third is the difference betwen B and A. The unique probability will
#                           be the number of elements in the set which is exclusive to A or 
#                           exclusive to B divided by the total number of peptides in the VIR 
#                           library pool. The exclusive set for A and the exclusive set for B
#                           could very well have different numbers, hence there needs to be an
#                           n^2 comparison.
#
#   Total probability   ->  This is the total peptides for each virus that has made it past the
#                           the bitscore of 80+. No comparisons are made to every other virus
#                           in the VIR lib peptide pool. Simply divide the number by the total
#                           count of peptides in VIR lib peptide pool.
#
 
def singleLoop(bitscoreDF, virlib_size, start):
    viruses = bitscoreDF.Species.drop_duplicates().values
    numVir = len(viruses)    
    unique = np.zeros(numVir)
    vA_set = set(
        bitscoreDF[bitscoreDF['Species'] == viruses[start]]['qseqid'].drop_duplicates().values
    )
    # print(f"Virus = {viruses[start]}, Number = {start}")
    for j in range(0, numVir):        
        vB_set = set(
            bitscoreDF[bitscoreDF['Species'] == viruses[j]]['qseqid'].drop_duplicates().values
        )
        exclusiveA = vA_set.difference(vB_set)        
        unique[j] = len(exclusiveA) / virlib_size        
    total = len(vA_set) / virlib_size
    uniqueDF = pd.DataFrame([unique], columns = viruses, index = [viruses[start]])
    totalDF = pd.DataFrame({'X1':[viruses[start]],'X2':[total]})
    return uniqueDF, totalDF


# ========================= group by taxon id ========================= #
# We group viral species by that taxon id. Then for each taxon id
# if there are multiple entries for their nt sequence, we only replace
# the nt sequence if it is longer than the one stored.
# After we've done this, write to file. This function is intended to 
# output file to be used with cd-hit-est so we can cluster viral species
# that are highly similar and reduce the number of comparisons needed
# when using AVARDA.

def groupByTaxonID(gb_dict,nt_output_prefix = None):
    taxon_id_dict = {}
    for accession_num in gb_dict.keys():
        taxon_id = gb_dict[accession_num].taxon_id
        organism = re.sub(' ','_',gb_dict[accession_num].organism)
        full_seq = gb_dict[accession_num].full_seq
        myKey = taxon_id + "|" + organism        
        if taxon_id in taxon_id_dict:
            if(len(taxon_id_dict[myKey]) < len(full_seq)):
                taxon_id_dict[myKey] = full_seq
        else:
            taxon_id_dict[myKey] = full_seq
    if(nt_output_prefix != None):
        nt_output = nt_output_prefix + ".fa"
        with open(nt_output, 'w') as fastaHandle:
            for key in taxon_id_dict.keys():
                fastaHeader = ">" + key
                s = taxon_id_dict[key]
                fastaHandle.write(fastaHeader + "\n" + s + "\n")
        fastaHandle.close()
    return taxon_id_dict


# ========================= read cdhit results ========================= #
# This function reads in the .clstr results from cd-hit-est. For each
# identified cluster, retreive the fasta description of the representative 
# sequence and return a list. The number of clusters the .clstr file 
# will be the length of the returned list.

def read_cdhit_clust(cdhit_path):    
    rep_seq_list = []    
    with open(cdhit_path, 'r') as cdhit_handle:
        for line in cdhit_handle:
            line = line.rstrip()
            if(line.startswith(">")):
                # New cluster
                pass
            elif (line.endswith("*")):
                repSeq = line.split(',')[-1]
                repSeq = re.sub('... \*','',re.sub(' >','',repSeq))
                rep_seq_list.append(repSeq)
    cdhit_handle.close()
    return rep_seq_list
            

# ========================= retrieving respresentative sequences ========================= #
# This function retrieves the representative sequence and writes a fasta file and a
# mapping file. The mapping file will be a tsv file that records taxon id to species
# name so we can retrieve it later.
# Note: Stupid makeblastdb does not allow continuous text to go above 50. Must add
#       space if longer.

def get_rep_seqs(rep_seq_list,taxon_id_dict,nt_output_prefix):
    taxonSeqKeys = taxon_id_dict.keys()
    nt_output = nt_output_prefix + ".fa"
    taxon_output = nt_output_prefix + ".tsv"
    with open(nt_output, 'w') as fastaHandle, open(taxon_output, 'w') as taxonHandle:
        taxonHandle.write( "sseqid" + "\t" + "Species" + "\n")
        for repSeq in rep_seq_list:
            # sanity check if statement - not really needed
            # if verified it works.
            if(repSeq in taxonSeqKeys):
                t = repSeq.split("|")[0] # taxon id
                v1 = repSeq.split("|")[1] # virus species
                v2 = re.sub('_',' ', v1)
                fastaHeader = ">" + t + " " + v2
                s = taxon_id_dict[repSeq]
                fastaHandle.write(fastaHeader + "\n" + s + "\n")
                taxonHandle.write( t + "\t" + v1 + "\n")
    fastaHandle.close()
    taxonHandle.close()    


# get_rep_seqs(testList, taxonStuff, "sequence.gb.testTaxonRepSeqs")
    
# ========================= constructing support files ========================= #
# tblastnOutput:    The output file from aligning our "nt_output" (above function)
#                   was used to build a viral database followed by the alignment
#                   through tblastn. We want to use the output.
# taxonMapFile:     The mapFile outputted from the above function.
#
# My example files:
# taxonMapFile -> sequence.gb.testTaxonRepSeqs.tsv
# tblastnOutput -> repSeq_tblastn_output
#
# If using WOOKIES_virlib_names_no91273.csv as the VIR library, then the virlib_size
# is 115752 peptides.
    
def createSupportTables(taxonMapFile, tblastnOutput, virlib_size, num_cores):    
    tblastn_80plusPath = tblastnOutput + "_bitscore80plus.tsv"
    tblastn_80minusPath = tblastnOutput + "_bitscore80minus.tsv"
    # read in taxon<->species map file
    taxonMap = dt.fread(
        taxonMapFile
    ).to_pandas()
    # read in tblastn output file and give columns their names
    repSeq_tblastn = dt.fread(
        tblastnOutput
    ).to_pandas().rename(
        columns = {
            "C0":"qseqid",
            "C1":"sseqid",
            "C2":"pident",
            "C3":"length",
            "C4":"mismatch",
            "C5":"gapopen",
            "C6":"qstart",
            "C7":"qend",
            "C8":"sstart",
            "C9":"send",
            "C10":"evalue",
            "C11":"bitscore"
        }
    )
    # filter out bitscore above 80 and merge to get species
    tblastn_80plus = repSeq_tblastn.loc[
        repSeq_tblastn['bitscore']>80
    ].merge(
        taxonMap, on = "sseqid"
    ).drop_duplicates().reset_index(drop = True)
    # filter out bitscore less than/equal to 80 and merge to get species
    # techincally we don't need this file but keeping it for reference
    tblastn_80minus = repSeq_tblastn.loc[
        repSeq_tblastn['bitscore']<=80
    ].merge(
        taxonMap, on = "sseqid"
    ).drop_duplicates().reset_index(drop = True)
    # print(tblastn_80plus)    
    # write out the files
    tblastn_80plus.to_csv(tblastn_80plusPath, index = False,  header = True,sep = '\t')
    tblastn_80minus.to_csv(tblastn_80minusPath, index = False, header = True, sep = '\t') 
    del(tblastn_80minus)   
    # Now we can calculate pident or bitscore table for 80plus only
    print(f"Making pep species tables...\n")
    pident = make_pep_species_table_v2('pident',tblastn_80plus)
    bitscore = make_pep_species_table_v2('bitscore',tblastn_80plus)
    bitscore.to_csv(tblastnOutput + "_bitscore80plus_bitscore.csv", index = True, header= True)
    pident.to_csv(tblastnOutput + "_bitscore80plus_pident.csv", index = True, header= True)
    # Next we calculate unique and total probability table for bitscore80plus
    # here we will use multiprocessing. Make sure num_core is set to
    # <= cpu_count() from multiprocessing.cpu_count
    assert(num_cores <= cpu_count()) 
    multiple_results = {}    
    pool = Pool(processes=num_cores)
    virus = tblastn_80plus.Species.drop_duplicates().values
    numVir = len(virus)
    print("Commence multiprocessing to calculate probabilities...\n")
    for i in range(numVir):
        single_result = pool.apply_async(
            singleLoop,
            [tblastn_80plus,virlib_size,i]
        )
        multiple_results[virus[i]] = single_result
        time.sleep(1)
    pool.close()
    pool.join()        
    print("Now the pool is closed and waiting for processes to finish for dev set...\n")
    # logging.info("Now the pool is closed and waiting for processes to finish for dev set...")
    u_collection = []
    t_collection = []
    for v in multiple_results.keys():
        single_result = multiple_results[v]
        u,t = single_result.get()
        u_collection.append(u)
        t_collection.append(t)
    totalDF = pd.concat(t_collection)    
    uniqueDF = pd.concat(u_collection)
    uniqueDF.to_csv(tblastnOutput + "_bitscore80plus_unique_probabilities_xr.csv")
    totalDF.to_csv(tblastnOutput + "_bitscore80plus_total_probabilities.csv", index = False)
    print(
        f"""
        Output files:
            *   Unique Probabilties:    {tblastnOutput}_bitscore80plus_unique_probabilities_xr.csv
            *   Total Probabilities:    {tblastnOutput}_bitscore80plus_total_probabilities.csv
            *   Pep2Pep (bitscore):     {tblastnOutput}_bitscore80plus_bitscore.csv
            *   Pep2Pep (pident):       {tblastnOutput}_bitscore80plus_pident.csv
        """
    )
    

# createSupportTables('sequence.gb.testTaxonRepSeqs.tsv','repSeq_tblastn_output',115752,4)
    
# ========================= count fasta size ========================= #
# count how many fasta sequences in a fasta file. Uses grep and wc -l
# because I am lazy.
    
def countFastaSize(fastafile):
    cmd1 = f"grep '>' {fastafile}"
    cmd2 = f"wc -l"
    output1 = sp.run(shlex.split(cmd1), capture_output=True)
    output2 = sp.run(shlex.split(cmd2), input=output1.stdout, capture_output=True)
    return int(output2.stdout.decode().strip())