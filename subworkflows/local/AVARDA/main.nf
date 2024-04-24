nextflow.enable.dsl=2
include { AVARDAOUT } from '../../../modules/local/avardaout/main.nf'


//
threshold_ch = Channel.value(params.threshold)
dict_ch = Channel.fromPath(params.dict_path)
total_ch = Channel.fromPath(params.total_path)
pairwse_ch = Channel.fromPath(params.pairwise_path)
blast_ch = Channel.fromPath(params.blast_path)
outpath_ch = Channel.value(params.)
outname_ch
workflow AVARDA{    
    take:
        // input 9 - avarda_names (path channel)
        virlib_ch
        // input 1 - case_path (path channel)
        edgeRhits_ch
    main:
        AVARDAOUT()
}

