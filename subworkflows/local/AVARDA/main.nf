nextflow.enable.dsl=2
include { AVARDAOUT } from '../../../modules/local/avardaout/main.nf'


if(params.run_AVARDA){
    threshold_ch = Channel.value(params.threshold)
    dict_ch = Channel.fromPath(params.dict_path)
    total_ch = Channel.fromPath(params.total_path)
    pairwse_ch = Channel.fromPath(params.pairwise_path)
    blast_ch = Channel.fromPath(params.blast_path)
    outpath_ch = Channel.value(params.out_path)
    outname_ch = Channel.value(params.out_name)
    cores_ch = Channel.value(params.max_cpus)
}

workflow AVARDA{    
    take:
        // input 9 - avarda_names (path channel)
        virlib_ch
        // input 1 - case_path (path channel)
        edgeRhits_ch
    main:
        AVARDAOUT(
            edgeRhits_ch,
            threshold_ch,
            dict_ch,
            total_ch,
            pairwse_ch,
            blast_ch,
            outpath_ch,
            outname_ch,
            virlib_ch,
            cores_ch
        )
}

