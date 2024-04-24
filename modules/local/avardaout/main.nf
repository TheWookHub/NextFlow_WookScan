process AVARDAOUT{
    publishDir "$params.out_path", mode: 'copy', overwrite: true
    input:
        path case_path
        val threshold
        path dict_path
        path total_path
        path pairwise_path
        path blast_path
        val out_path
        val out_name
        path avarda_names
    
    script:
    """
    echo 
        case_path:          $case_path
        avarda_names:       $avarda_names
    """
        // '''
        //     Rscript AVARDA.R \
        //     --case_path /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/subworkflows/local/AVARDA/data/example_input/AVARDA_test_data.tsv.gz \
        //     --threshold 1 \
        //     --dict_path /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/subworkflows/local/AVARDA/data/dict_path/my_df.csv.gz \
        //     --total_path /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/subworkflows/local/AVARDA/data/total_path/total_probability_xr2.csv \
        //     --pairwise_path /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/subworkflows/local/AVARDA/data/pairwise_path/unique_probabilities3.csv \
        //     --blast_path /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/subworkflows/local/AVARDA/data/blast_path/VirScan_filtered_virus_blast_new.csv.gz \
        //     --out_path /home/preston/PhIPSeq-Pipelines/AVARDA_Runs/TestWookflow_AVARDA \
        //     --out_name WATERMELON  \
        //     --avarda_names /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/subworkflows/local/AVARDA/data/avarda_names/avarda_names.csv.gz
        // '''
}
