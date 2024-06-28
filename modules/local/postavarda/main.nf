process POST_AVARDA_PROCESSING{
    publishDir "$params.out_path", mode: 'copy', overwrite: true
    container = 'quay.io/hdc-workflows/phippery:1.2.0'
    debug true
    input:
        path sample_table
        val file_prefix                
    output:
        path "*.csv", emit: postavarda_outs        
    script:    
    """
    avardaSampleID.py \
        -i $sample_table \
        -p $file_prefix
    """
}
