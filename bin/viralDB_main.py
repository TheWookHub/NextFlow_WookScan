#!/usr/bin/python


#==========================="Python library imports"===========================#
import sys, argparse, os, shlex, psutil, math
import subprocess as sp
import viralDB as vDB
import viralDB_classes as vDBc
from multiprocessing import cpu_count

#==========================="Global Variables"===========================#

WORKDIR = os.getcwd()
MAX_CORES = cpu_count()
NUM_CORES = 4
FILE_PREFIX = "viraldb_output"
INTERNALDATA = "ViralDB_Data"
VIRLIB = os.path.join(INTERNALDATA,"VIR3_Peptide_noDup.fa")
NUM_ID = 115752
CDCLUSTCUTOFF = 0.8
CDMEM = 16000 #in Megabytes
MAX_SYSMEM = psutil.virtual_memory().total / math.pow(1024,2) # Max system memory in Megabytes

#==========================="Main function"===========================#
if __name__ == '__main__':

    program_function = """
    *****
    ViralDB Tool:
        This tool takes in a genbank file (.gb) consisting of virus entries and a 
        viralhostdb (.tsv) to create the supporting files for AVARDA.

        *   The genbank file can be obtained at:
                https://www.ncbi.nlm.nih.gov/genome/viruses/
            via selecting "RefSeq and Neighbour nucleotide records" options.
            However this resource seems to be discontinuing. In that case it also works
            if you paste:
        
            "
                Viruses[Organism] NOT cellular organisms[ORGN] NOT wgs[PROP] NOT AC_000001:AC_999999[pacc]
                NOT gbdiv syn[prop] AND (srcdb_refseq[PROP] OR nuccore genome samespecies[Filter]) 
                AND ("vhost human"[Filter]) 
            "

            As search criteria in genbank.

        *   The .tsv file from viralhostdb can be obtained at:
                https://www.genome.jp/virushostdb/
            and select Human Viruses.

        *   This tool currently only looks for entries that have human (homo sapiens)
            as host.
        
    *****
    """
    parser = argparse.ArgumentParser()

    # Required Arguments
    parser.add_argument('-g', '--genBankFile', help = 'Filename of the genebank (.gb) file.',  required = True)
    parser.add_argument('-v','--vhdb', help="Viral Host Database file (.tsv)", required=True)    
    
    # Optional Arguments
    parser.add_argument('-cpu','--processors', help = f"Number of processes for parallel computation. [Default = {NUM_CORES}]", type = int)    
    parser.add_argument('-o','--out', help = f"Output path. [Default:{WORKDIR}]", type = str)   
    parser.add_argument('-p','--prefix',help = f"File prefix for outputfiles. [Default: {FILE_PREFIX}]", type = str)
    parser.add_argument('-l','--virlib',help=f"Viral peptide library (.fa) used in PhIPSeq. [Default: {VIRLIB}]")
    parser.add_argument('-cc','--cluster_cutoff',help=f" [Default: {CDCLUSTCUTOFF}]", type = float)
    parser.add_argument('-cm','--cluster_mem', help=f"[Default: {CDMEM}]", type = int)

    if(len(sys.argv) < 2):
        print(program_function)
        parser.print_help()
        sys.exit()
    args=parser.parse_args()
    
    # check if the input files are valid    
    if(not os.path.isfile(args.genBankFile)):
        raise vDBc.InputFileError(args.genBankFile,  "Can't find file: ")
    if(not os.path.isfile(args.vhdb)):
        raise vDBc.InputFileError(args.vhdb,  "Can't find file: ")    
    if(args.out):
        if(not os.path.exists(args.out)):
            print(f"{args.out} not found! Creating directory...")
            os.makedirs(args.out)        
        WORKDIR = args.out
    if(args.virlib):
        if(not os.path.isfile(args.virlib)):
            raise vDBc.InputFileError(args.virlib,  "Can't find file: ")
        VIRLIB = args.virlib
        NUM_ID = vDB.countFastaSize(VIRLIB)
    if(args.prefix):
        FILE_PREFIX = args.prefix
    if(args.processors > MAX_CORES):
        print(f"This system only has {MAX_CORES} cores. Given {args.processors}. Revert to {MAX_CORES - 1}.")
        NUM_CORES = MAX_CORES - 1
    else:
        NUM_CORES = args.processors
    if(args.cluster_cutoff):
        if(args.cluster_cutoff < 0 | args.cluster_cutoff > 100):
            raise vDBc.InputFileError(args.cluster_cutoff, "Value out of bounds (between 0.0 and 100.0). Given: ")
        else:
            CDCLUSTCUTOFF = args.cluster_cutoff
    if(args.cluster_mem):
        if(args.cluster_mem > MAX_SYSMEM):
            print(f"Warning: cluster_mem {args.cluster_mem} MB is > System's available memory: {MAX_SYSMEM} MB.")
            CDMEM = MAX_SYSMEM * 0.9
            print(f"Reducing the memory request to {CDMEM} MB.")
        else:
            CDMEM = args.cluster_mem


    print(
        f"# ---------------------------------- User given parameters: ---------------------------------- #"
    )
    print(
        f"""        
            Genbank File:               {args.genBankFile}
            Viral Host Database File:   {args.vhdb}
            VirLib Used:                {VIRLIB}
            No. Peptide IDs             {NUM_ID}
            No. CPU:                    {NUM_CORES}
            Requested Memory:           {CDMEM}
            Clustering Cutoff:          {CDCLUSTCUTOFF}
            Output Dir:                 {WORKDIR}
            File Prefix:                {FILE_PREFIX}
        """
    )
    
    print(
        f"# ------- Step 1: Reading in files ------- #\n"
    )
    vhdb = vDB.getVHDB_humans(args.vhdb)
    myGB = vDB.getGenBankInfo(args.genBankFile, vhdb)

    print(
        f"# ------- Step 2: Grouping by taxon id ------- #\n"
    )
    outputname = os.path.join(WORKDIR,FILE_PREFIX)
    # writes viraldb_output.fa (default file name) ready for cdhit 
    taxon_id_dict = vDB.groupByTaxonID(myGB,outputname) 

    print(
        f"# ------- Step 3: CD-HIT-EST to retrieve representative taxonmy & sequences ------- #\n"
    )    
    # reads in viraldb_output.fa (default file name)
    # outputs viraldb_output_90Clust.clstr
    tool = "cd-hit-est"
    params1 = f" -i {outputname}.fa" 
    params2 = f" -o {outputname}_Clust"
    params3 = f" -c {CDCLUSTCUTOFF} -M {CDMEM} -T {NUM_CORES} -d 0"
    sp.run(shlex.split(tool + params1 + params2 + params3))
    rep_seq_list = vDB.read_cdhit_clust(outputname +'_Clust.clstr')    
    # viraldb_output_rep
    rep_output = f"{outputname}_rep"
    rep_output_fname = f"{rep_output}.tsv"
    vDB.get_rep_seqs(rep_seq_list,taxon_id_dict,rep_output)
    print(f"")
    print(
        f"# ------- Step 4: building viral db, tblastn and blastp (params set according to AVARDA paper) ------- #\n"
    )
    # These args are hardcoded because it is following the settings from AVARDA paper
    # viraldb_output_rep.fa
    
    tool = "makeblastdb" 
    params = f" -in {rep_output}.fa -parse_seqids -dbtype 'nucl'"
    print(f"cmd:    {tool + params}")
    sp.run(shlex.split(tool + params))   
    print(f"{tool}: completed\n")
    
    tool = 'tblastn'
    tblastn_fname = f"{rep_output}_tblastn_output"
    params1 = f" -query {VIRLIB} -db {rep_output}.fa -out {rep_output}_tblastn_output"
    params2 = f" -word_size 7 -outfmt 6 -seg no -soft_masking false -max_hsps 1 -max_target_seqs 100000" 
    params3 = f" -num_threads {NUM_CORES}"
    print(f"cmd:    {tool + params1 + params2 + params3}")
    sp.run(shlex.split(tool + params1 + params2 + params3))
    print(f"{tool}: completed\n")
    
    # attempt to split the fasta file into how many CPU there are
    
    VIRLIB_SPLITS = vDB.createFastaSplits(NUM_CORES,VIRLIB)
    blastp_fname = vDB.run_blastp(VIRLIB,VIRLIB_SPLITS,NUM_CORES,outputname)
    # later on for improvement - we can make this into a loop
    # split the VIRLIB file into N number of allowed cores
    # and then run in parallel.
    # Note:     -query would be the not split file and
    #           -subject would change when parallel runing. One
    #           for each split.

    # tool = 'blastp'
    # blastp_fname = f"{outputname}_blastp_output"
    # params1 = f" -query {VIRLIB} -subject {VIRLIB}"
    # params2 = f" -out {blastp_fname}"
    # params3 = f" -outfmt 6 -evalue 100"
    # print(f"cmd:    {tool + params1 + params2 + params3}")
    # sp.run(shlex.split(tool + params1 + params2 + params3))
    # print(f"{tool}: completed\n")

    print(
        f"# ------- Step 5: Creating AVARDA support files ------- #\n"
    )
    vDB.createSupportTables(rep_output_fname, tblastn_fname, blastp_fname, NUM_ID, NUM_CORES)
    # rubbish = VIRLIB_SPLITS
    # vDB.cleanUp(rubbish)
    print(
        f"# ------- Done! ------- #"
    )