process AVARDAOUT{
    publishDir "$params.out_path", mode: 'copy', overwrite: true
    label 'process_long'
    container = 'docker.io/pdawgzgg/avarda_r_env:0.1'
    debug true
    input:
        path case_path
        val threshold
        path dict_path
        path total_path
        path pairwise_path
        path blast_path
        val out_name
        path avarda_names
        val cores
    output:
        path "*.csv", emit: avarda_outs        
    script:    
    """
    AVARDA.R \
        --case_path $case_path \
        --threshold $threshold \
        --dict_path $dict_path \
        --total_path $total_path \
        --pairwise_path $pairwise_path \
        --blast_path $blast_path \
        --out_name $out_name \
        --avarda_names $avarda_names \
        --cores $cores
    """
}
