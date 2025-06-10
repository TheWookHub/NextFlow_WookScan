// process FASTP_OUT is to take in the fastq 

// probably don't need this script as I've already put it as a process
// in the fastp workflow.

// process FASTP_OUT{
//     publishDir "$params.results/pickle_data/", mode: 'copy', overwrite: true
//     input:
//         val upep_string
//         val phipdata_name
//     output:
//         path "*_virlib_names.csv", emit: virlib
//         path "PhipperyEdgeRHITS_AVARDA_Input.csv", emit: edgeRhits
//     script:
//         """
//         phippery_to_avarda_00.py \
//         -i $phipdata_name \
//         -u_pep_id $upep_string
//         """
// }
