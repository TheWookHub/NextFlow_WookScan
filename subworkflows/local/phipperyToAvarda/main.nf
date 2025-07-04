/*
*   This workflow is to prepare Phippery output so it can be used as input file
*   for AVARDA
*/

nextflow.enable.dsl=2

/* 
*   Create:
*       1)  virlib table 
*       2)  phippery edgeR avarda input read files
*/
// include { PHIPOUTPUT } from '../../../modules/local/phipout/main.nf'

/* 
*   If needed we generate the new viral database for avarda creating a series
*   of files. These include:
*       1)  New network file based on blastp peptide to peptide e-value score < 100. These
*           will be collected as edges in a csv file. 
*       2)  Generate bitscore80plus files:
*           -   total probabilities (species to species probability)
*           -   unique probabilities (species probability)
*           -   bitscore80plus bitscore table (peptide id to species bitscore)
*           -   bitscore80plus pident table (peptide id to species pident) <--- likely not needed but optional
*
*   This component requires the user to have downloaded a .gb file from genbank and a .csv (or .tsv) file
*   from ViralHostDB to feed in to this nexflow component.
*
*   Note to self: Need to figure out later how to co-ordinate the options to run this nexflow pipeline.
*/
// include { VIRALDB } from '../../../modules/local/phipout/main.nf'


process PHIPOUTPUT{
    publishDir "$params.results/pickle_data/", mode: 'copy', overwrite: true
    input:
        val upep_string
        val phipdata_name
    output:
        path "*_virlib_names.csv", emit: virlib
        path "PhipperyEdgeRHITS_AVARDA_Input.csv", emit: edgeRhits
    script:
        """
        phippery_to_avarda_00.py \
        -i $phipdata_name \
        -u_pep_id $upep_string
        """
}

process PHIPOUTPUT_EXTRACT{
    publishDir "$params.results/pickle_data/", mode: 'copy', overwrite: true
    input:        
        val phipdata_name
    output:
        path "*.csv"
    script:
        """
        simple_phippery_process.py \
        -d $phipdata_name

        """
}


if(params.run_phippery){
    upep_prefix_ch = Channel.value(params.user_pep_id)
}

workflow PHIPPERYTOAVARDA{
    take:
        data_phip_ch

    main:
        PHIPOUTPUT(upep_prefix_ch,data_phip_ch)
        PHIPOUTPUT_EXTRACT(data_phip_ch)
    emit:
        virlib = PHIPOUTPUT.out.virlib
        edgeRhits = PHIPOUTPUT.out.edgeRhits    
}
