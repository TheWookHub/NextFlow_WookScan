/*
*   This work flow is to do some post AVARDA Processing
*   such as re-linking the fastq files back to AVARDA output ids (it can be confusing)
#   and any other tools (potential) that help interpret AVARDA outputs
*/

nextflow.enable.dsl=2

/* 
*   Creates:
*       1)  an AVARDA 'name' to fastq file path csv so we know which result
*           is relevant to which input fastq file
*/
include { POST_AVARDA_PROCESSING } from '../../../modules/local/postavarda/main.nf'

// Grab the required param inputs
sample_table = Channel.fromPath(params.sample_table)
file_prefix = Channel.value(params.out_name)

workflow POSTAVARDA_WORKFLOW{    
    POST_AVARDA_PROCESSING(sample_table,file_prefix)
}


