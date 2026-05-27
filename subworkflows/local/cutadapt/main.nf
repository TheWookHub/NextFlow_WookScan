/*
*   WIP:
*   This nf script is to work with cutadapt. The goal is the trim 5' ends of the reads to remove the 5' adapters with 
*   variable lengths. These specific adapters can be between 3 to 6 in length.
*/

nextflow.enable.dsl=2

// process CUTADAPT_OUT is to take in the fastq and then run cutadapt for trimming
// we will rename these files with the "trimmed" suffix. 
process CUTADAPT_OUT{    
    publishDir "${params.results}/trimmed_fastq/", mode: 'copy', overwrite: true
    // container = 'docker.io/pdawgzgg/avarda_r_env:0.1'
    input:
        tuple val(tech_id),val(basename_r1),path(filename_r1),val(basename_r2),path(filename_r2)
    output:
        path("*_trimmed.fastq.gz"), emit: trimmedFqName        

    // need to change this cutadapt file input line because now there are two files instead of one
    // need to figure out how the input sample table file should look like
    script:
        """
        cutadapt \
        -j ${task.cpus} \
        -g "file:${params.r1_adapters}" \
        -G "file:${params.r2_adapters}" \
        -o "${basename_r1}_trimmed.fastq.gz" \
        -p "${basename_r2}_trimmed.fastq.gz" \
        ${filename_r1} ${filename_r2}
        """
}

process UNRAVEL_TRIMMED_NAMES{
    publishDir "$params.results/trimmed_fastq/", mode: 'copy', overwrite: true
    tag "PE unravel trimmed names"
    input:
        val trimmed_fastq_list
    output:
        path "trimmed_names_collection.tsv", emit: trimmed_names_collection
    script:
        """
        echo ${trimmed_fastq_list} | \
        grep -o '/[^ ]*fastq.gz' | \
        awk '/_R1_/ {r1=\$0} /_R2_/ {print r1 "\t" \$0}' > trimmed_names_collection.tsv

        """
}

// process UPDATE_SAMPLE_TABLE is to update the sample table with the filtered fastq file names
process UPDATE_SAMPLE_TABLE{
    tag "PE sample table update"
    publishDir "$params.results/trimmed_fastq/", mode: 'copy', overwrite: true
    input:
        path trimmed_fastq_csv
        path sample_table        
    output:
        path "trimmed_sample_table.csv", emit: trimmed_table
    script:
        """
        update_sample_table_huscan.py \
        -s ${sample_table} \
        -t ${trimmed_fastq_csv} \
        -r ${params.results}/trimmed_fastq \
        -o "trimmed_sample_table.csv"
        """
}

workflow CUTADAPT_WORKFLOW{
    take:
        sample_ch
    main:
        // Take original sample table
        // got to extract the first column that contains all
        // the fastq file paths and then run fastp.
        //sample_ch = Channel.fromPath(params.sample_table)
        sample_ch
            .splitCsv(header:true)
            .map{ row -> 
                    tuple(
                        row.technical_replicate_id,                         
                        (row.fastq_filepath =~ /.+\/(.+)\.fastq\.gz/)[0][1], // extract the base name of r1 fastq file
                        file(row.fastq_filepath), // fastq file path for r1
                        (row.fastq_filepath2 =~ /.+\/(.+)\.fastq\.gz/)[0][1], // extract the base name of r2 fastq file
                        file(row.fastq_filepath2) // fastq file path for r2
                    )                    
                }
            .set { sample_table_ch }        
        
        // Run cutadapt to remove 5' adapters in paired end fastq files.
        // then record the trimmed files and prep it to be pass to fastp.
        // NOTE: NEED TO CHANGE SO FASTP CAN TAKE IN PAIRED END FILES AT THIS STAGE
        CUTADAPT_OUT(sample_table_ch)        
        UNRAVEL_TRIMMED_NAMES(CUTADAPT_OUT.out.trimmedFqName.toList())        
        UNRAVEL_TRIMMED_NAMES.out.trimmed_names_collection.set{collection_tsv_ch}
        UPDATE_SAMPLE_TABLE(collection_tsv_ch,sample_ch)        
    emit:       
        sample_info = UPDATE_SAMPLE_TABLE.out.trimmed_table
}

