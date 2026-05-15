/*
*   WIP:
*   This nf script is to work with cutadapt. The goal is the trim 5' ends of the reads to remove the 5' adapters with 
*   variable lengths. These specific adapters can be between 3 to 6 in length.
*/

nextflow.enable.dsl=2

// process CUTADAPT_OUT is to take in the fastq and then run cutadapt for trimming
// we will rename these files with the "trimmed" suffix. 
process CUTADAPT_OUT{
    publishDir "$params.results/filtered_fastq/", mode: 'copy', overwrite: true
    //container = 'docker.io/pdawgzgg/avarda_r_env:0.1'
    input:
        tuple val(tech_id),val(basename),val(filename), path(file_path)
    output:
        path("*_filtered.fastq.gz"), emit: filteredFqName
        path ("*_filtered.html")

    // need to change this cutadapt file input line because now there are two files instead of one
    // need to figure out how the input sample table file should look like
    script:
        """        
        fastp -t 0 \
        -i "${file_path}" \
        -z 9 \
        -o "${basename[0][1]}_filtered.fastq.gz" \
        -R "${basename[0][1]}" \
        -j "${basename[0][1]}_filtered.json" \
        -h "${basename[0][1]}_filtered.html"

        cutadapt \
        -j "${task.cpus}" \
        -g "file:${params.r1_adapters}" \
        -G "file:${params.r2_adapters}" \
        -o "${basename[0][1]}_trimmed.fastq.gz" \
        -p "${basename[0][1]}_trimmed.fastq.gz" \
        HuScanPhageExp-V5-3_S1_L001_R1_001.fastq.gz HuScanPhageExp-V5-3_S1_L001_R2_001.fastq.gz
        """
}

process UNRAVEL_FILTERED_NAMES{
    publishDir "$params.results/filtered_fastq/", mode: 'copy', overwrite: true
    input:
        val filtered_fastq_list
    output:
        path "filtered_names_collection.txt", emit: filtered_names_collection
    script:
        """
        for item in ${filtered_fastq_list.join(' ')}; 
        do            
            echo "\$item" >> filtered_names_collection.txt
        done

        """
}

// process UPDATE_SAMPLE_TABLE is to update the sample table with the filtered fastq file names
process UPDATE_SAMPLE_TABLE{
    // publishDir "$params.results/filtered_fastq/", mode: 'copy', overwrite: true
    input:
        path filtered_fastq_list
        path sample_table        
    output:
        path "filtered_sample_table.csv", emit: filtered_table
    script:
        """
        update_sample_table.py \
        -s ${sample_table} \
        -t ${filtered_fastq_list} \
        -r ${params.results}/filtered_fastq \
        -o "filtered_sample_table.csv"
        """
}

workflow CUTADAPT_WORKFLOW{
    main:
        // Take original sample table
        // got to extract the first column that contains all
        // the fastq file paths and then run fastp.
        sample_ch = Channel.fromPath(params.sample_table)
        sample_ch
            .splitCsv(header:true)
            .map{ row -> 
                    tuple(
                        row.technical_replicate_id, 
                        row.fastq_filepath =~ /.+\/(.+)\.fastq\.gz/, // extract the base name of the fastq file
                        row.fastq_filepath, // fastq file path
                        row.fastq_filepath2 =~ /.+\/(.+)\.fastq\.gz/, // extract the base name of the fastq file
                        row.fastq_filepath2, // fastq file path
                        file("$params.reads_prefix/${row.fastq_filepath}")
                    ) 
                }
            .set { sample_table_ch }        
        
        // Run fastp for quality check and then collect the new filtered fastq
        // file names. After that alter the original sample table so the fastq file paths
        // are pointing to the new filtered fastq files.        
        // CUTADAPT_OUT(sample_table_ch)         
        // UNRAVEL_FILTERED_NAMES(FASTP_OUT.out.filteredFqName.toList())
        // UNRAVEL_FILTERED_NAMES.out.filtered_names_collection.set{collection_ch}
        // UPDATE_SAMPLE_TABLE(collection_ch,sample_ch)        
    emit:
        sample_info = sample_table_ch
    //     sample_info = UPDATE_SAMPLE_TABLE.out.filtered_table
}

