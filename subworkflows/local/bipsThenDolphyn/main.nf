#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/bipsdolmethodtwo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/bipsdolmethodtwo
    Website: https://nf-co.re/bipsdolmethodtwo
    Slack  : https://nfcore.slack.com/channels/bipsdolmethodtwo
----------------------------------------------------------------------------------------
*/

/*
*   This workflow is to process four different modes of BIPS and Dolphyn.
*   They are BIPS-Related: bips_then_dolphyn, bips_only.
*   Dolphyn-related: dolphyn_standalone, dolphyn_oligo_only.
*   
*/

nextflow.enable.dsl=2




// --- Processes ---

process RUN_BIPS_INITIAL {
    // The tag now represents a batch run, so sample_id is not directly applicable.
    // We can use a fixed tag or one based on the mode.
    tag "BIPS Viral batch run for mode ${params.mode}"

    publishDir "${params.outdir}/bips_batch_run_outputs", // Note: Fixed name for batch mode
        mode: 'copy',
        overwrite: true    

    input:
    path file_list
    path bips_code_dir

    output:
    // These paths are now relative to the Nextflow work directory after being copied from BIPS_run
    // They are a single set of results for the entire batch.
    path(params.bips_oligos_sequence_csv_name), emit: oligos_csv
    path(params.bips_barcoded_nuc_file_csv_name), emit: barcoded_csv
    path "bips_full_output" 

    script:
    def bips_staged_name = bips_code_dir.baseName // e.g., "BuildPhIPSeqLibrary"
    def oligos_out_name = params.bips_oligos_sequence_csv_name
    def barcoded_out_name = params.bips_barcoded_nuc_file_csv_name
    def full_output_container_name = "bips_full_output" 

    """
    echo "RUN_BIPS_INITIAL for a batch of files"
    echo "Staged BIPS code directory: ${bips_staged_name}"

    # --- Clean BIPS I/O directories using Shell commands ---
    BIPS_INPUT_DIR="${bips_staged_name}/Data/Input"
    BIPS_OUTPUT_DIR="${bips_staged_name}/Data/Output"

    echo "Cleaning \$BIPS_INPUT_DIR (keeping README.md)..."
    mkdir -p "\$BIPS_INPUT_DIR" # Ensure it exists
    find "\$BIPS_INPUT_DIR" -mindepth 1 -type f ! -name 'README.md' -delete
    find "\$BIPS_INPUT_DIR" -mindepth 1 -type d -empty -delete # Remove empty subdirs if any

    echo "Cleaning \$BIPS_OUTPUT_DIR (keeping README.md)..."
    mkdir -p "\$BIPS_OUTPUT_DIR" # Ensure it exists
    find "\$BIPS_OUTPUT_DIR" -mindepth 1 -type f ! -name 'README.md' -delete
    find "\$BIPS_OUTPUT_DIR" -mindepth 1 -type d -empty -delete
    echo "BIPS I/O directories cleaned."
    # --- End Cleaning Step ---


    # --- Prepare BIPS input by copying ALL files from the list ---

    echo "Copying multiple input files to BIPS input directory..."
    for file in ${file_list}; do
        # For each file in the list, copy it to the BIPS input directory
        cp "\$file" "\$BIPS_INPUT_DIR/"
        echo "  Copied \$file"
    done
    echo "All input files copied."


    # --- Run BIPS ---
    echo "Changing to BIPS directory: ${bips_staged_name}"
    cd "${bips_staged_name}"

    echo "Listing Data/Input/:"
    ls -l Data/Input/

    echo "Executing BIPS: python3 main.py"
    python3 main.py > bips_stdout.log 2> bips_stderr.log
    BIPS_EXIT_CODE=\$?

    cd .. # Return to original work directory

    echo "--- BIPS STDOUT ---"
    cat "${bips_staged_name}/bips_stdout.log"
    echo "--- BIPS STDERR ---"
    cat "${bips_staged_name}/bips_stderr.log" >&2

    # --- Check Exit Code and File Existence ---
    if [ "\$BIPS_EXIT_CODE" -ne 0 ]; then
        echo "ERROR: BIPS main.py failed with exit code \$BIPS_EXIT_CODE" >&2
        exit \$BIPS_EXIT_CODE
    fi

    # === ADD THE COPY COMMANDS BACK IN ===
    # Copy results from BIPS's output dir to the root of the work dir
    # This makes the files available for the 'output:' block declarations.

    BIPS_RAW_OUTPUT_DIR="${bips_code_dir}/Data/Output"

    # Check if the BIPS output directory itself was created
    if [ ! -d "\$BIPS_RAW_OUTPUT_DIR" ]; then
        echo "ERROR: BIPS did not create the output directory: \$BIPS_RAW_OUTPUT_DIR" >&2
        exit 1
    fi

    mkdir -p "${full_output_container_name}"


    # 3. Copy ALL of BIPS's results into this new directory
    echo "Copying all BIPS results to '${full_output_container_name}' for publishing..."
    if [ -n "\$(ls -A \$BIPS_RAW_OUTPUT_DIR)" ]; then
        cp -rL "\$BIPS_RAW_OUTPUT_DIR"/* "${full_output_container_name}/"
    else
        echo "WARNING: BIPS ran successfully but produced no output files in \$BIPS_RAW_OUTPUT_DIR"
    fi


    # 4. Copy the two key files needed by downstream channels to the work dir root
    if [ -f "${full_output_container_name}/${oligos_out_name}" ]; then
        cp "${full_output_container_name}/${oligos_out_name}" .
    else
        echo "ERROR: BIPS output ${oligos_out_name} not found in results!" >&2; exit 1
    fi
    if [ -f "${full_output_container_name}/${barcoded_out_name}" ]; then
        cp "${full_output_container_name}/${barcoded_out_name}" .
    else
        echo "ERROR: BIPS output ${barcoded_out_name} not found in results!" >&2; exit 1
    fi
    """
}

process BIPS_CSV_TO_FASTA{
    tag "${sample_id}"
    publishDir "${params.outdir}/dolphyn_prep", mode: 'copy', pattern: "*.fasta"

    input:
    tuple val(sample_id), path(oligos_csv)
    
    output:
    tuple val(sample_id), path("${sample_id}.oligos_for_dolphyn.fasta"), emit: oligos_fasta
    
    when:
    (params.mode == 'bips_then_dolphyn') ||
    (params.mode == 'dolphyn_oligo_only' && params.input_oligos_csv_dir)
    
    script:
    """
    ${params.helper_script} bips_csv_to_fasta \
        --bips_oligos_csv ${oligos_csv} \
        --output_fasta ${sample_id}.oligos_for_dolphyn.fasta
    """
}

process RUN_DOLPHYN_PREDICTION {
    conda "$baseDir/bips_environment.yml"
    
    tag "$sample_id"
    publishDir "${params.outdir}/dolphyn_prediction", mode: 'copy', pattern: "*.json"

    input:
    tuple val(sample_id), path(oligos_fasta) // from BIPS_CSV_TO_FASTA

    output:
    tuple val(sample_id), path("${sample_id}.dolphyn_raw.json"), emit: dolphyn_json

    script:
    // params.dolphyn_package_path should be "${projectDir}/vendor/Dolphyn"
    // This is the directory containing the 'dolphyn' Python package.
    def dolphyn_path = params.dolphyn_package_path

    // export PYTHONPATH="${params.dolphyn_package_path}:\${PYTHONPATH}"

    """
    # Safely set or prepend to PYTHONPATH
    if [ -z "\${PYTHONPATH:-}" ]; then  # Check if PYTHONPATH is unset or empty
      export PYTHONPATH="${dolphyn_path}"
    else
      export PYTHONPATH="${dolphyn_path}:\${PYTHONPATH}"
    fi
    echo "PYTHONPATH set to: \$PYTHONPATH" # For debugging
    
    ${params.helper_script} run_dolphyn \\
        --input_fasta ${oligos_fasta} \\
        --output_json ${sample_id}.dolphyn_raw.json \\
        --dolphyn_training_data_dir ${params.dolphyn_training_data_dir}
    """
}

process DOLPHYN_JSON_TO_CSV {
    tag "$sample_id"
    publishDir "${params.outdir}/dolphyn_conversion", mode: 'copy', pattern: "*.csv"

    input:
    tuple val(sample_id), path(dolphyn_json) // from RUN_DOLPHYN_PREDICTION

    output:
    tuple val(sample_id), path("${sample_id}.epitopes_from_dolphyn.csv"), emit: epitope_csv

    script:
    """
    ${params.helper_script} dolphyn_json_to_csv \\
        --dolphyn_json ${dolphyn_json} \\
        --epitope_csv_output ${sample_id}.epitopes_from_dolphyn.csv
    """
}

process SELECT_EPITOPE_POSITIVE_BIPS_OLIGOS {
    tag "$sample_id"
    publishDir "${params.outdir}/matching_step", mode: 'copy', pattern: "*.csv"

    input:
    // original_bips_oligos_csv comes from RUN_BIPS_INITIAL
    // dolphyn_epitopes_csv comes from DOLPHYN_JSON_TO_CSV
    tuple val(sample_id), path(original_bips_oligos_csv), path(dolphyn_epitopes_csv)

    output:
    tuple val(sample_id), path("${sample_id}.selected_bips_oligos.csv"), emit: selected_bips_oligos_csv

    script:
    """
    ${params.helper_script} select_epitope_positive_bips_oligos \\
        --original_bips_oligos_csv ${original_bips_oligos_csv} \\
        --dolphyn_epitopes_csv ${dolphyn_epitopes_csv} \\
        --selected_bips_oligos_output_csv ${sample_id}.selected_bips_oligos.csv
    """
}

process FILTER_BIPS_BARCODES {
    tag "$sample_id"
    publishDir "${params.outdir}/final_filtered_barcodes", mode: 'copy', pattern: "*.csv"

    input:
    // full_bips_barcoded_csv comes from RUN_BIPS_INITIAL
    // selected_bips_oligos_csv comes from SELECT_EPITOPE_POSITIVE_BIPS_OLIGOS
    tuple val(sample_id), path(full_bips_barcoded_csv), path(selected_bips_oligos_csv)

    output:
    tuple val(sample_id), path("${sample_id}.final_epitope_barcoded_oligos.csv"), emit: filtered_barcoded_csv

    script:
    """
    ${params.helper_script} filter_bips_barcodes \\
        --full_barcoded_csv ${full_bips_barcoded_csv} \\
        --selected_bips_oligos_csv ${selected_bips_oligos_csv} \\
        --filtered_barcoded_output_csv ${sample_id}.final_epitope_barcoded_oligos.csv
    """
}

process RUN_DOLPHYN_STANDALONE_PREP {
    conda "$baseDir/bips_environment.yml"


    tag "$sample_id (Dolphyn Standalone)"
    publishDir "${params.outdir}/dolphyn_standalone_prediction", mode: 'copy', pattern: "*.json"

    input:
    tuple val(sample_id), path(protein_fasta_file)
    // path dolphyn_training_data_dir from file(params.dolphyn_training_data_dir) // If still needed by action

    output:
    tuple val(sample_id), path("${sample_id}.dolphyn_standalone.json"), emit: dolphyn_json

    script:
    def dolphyn_path = params.dolphyn_package_path

    """
    # Safely set or prepend to PYTHONPATH
    if [ -z "\${PYTHONPATH:-}" ]; then  # Check if PYTHONPATH is unset or empty
      export PYTHONPATH="${dolphyn_path}"
    else
      export PYTHONPATH="${dolphyn_path}:\${PYTHONPATH}"
    fi
    echo "PYTHONPATH set to: \$PYTHONPATH" # For debugging


    ${params.helper_script} run_dolphyn_standalone \\
        --input_protein_fasta ${protein_fasta_file} \\
        --output_json ${sample_id}.dolphyn_standalone.json \\
        --dolphyn_training_data_dir ${params.dolphyn_training_data_dir} 
    """
}



// --- Workflow ---
workflow BIPS_THEN_DOLPHYN {
    take:
        ch_input_for_bips  // ch_viral_seqs_for_bips		// For modes: bips_then_dolphyn, bips_only
        ch_input_for_dolphyn_standalone // ch_viral_seqs_for_dolphyn_standalone		// For mode: dolphyn_standalone
        ch_oligos_fasta_for_dolphyn_only		// For mode: dolphyn_oligo_only (oligo FASTA)
        ch_oligos_csv_for_dolphyn_only		// For mode: dolphyn_oligo_only (oligo CSV)


        //ch_input_files
        bips_root_path_obj // Path object for BIPS vendored code
        helper_script_path_obj // Path object for your Python helper

    main:
			// == Initialize Core Channels ==
			// These will be populated based on the pipeline mode.
			ch_bips_oligos_csv_result = Channel.empty()
			ch_bips_barcoded_csv_result = Channel.empty()

			ch_fasta_for_core_dolphyn = Channel.empty() // This is for oligo-based Dolphyn prediction
			ch_dolphyn_json_from_core_prediction = Channel.empty()
			ch_dolphyn_epitope_csv_result = Channel.empty()
			ch_selected_bips_oligos_result = Channel.empty()
			ch_final_filtered_barcodes_result = Channel.empty()
			ch_dolphyn_standalone_json_result = Channel.empty()


			// ========================
			//      BIPS Execution
			// ========================
			if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only") {
					RUN_BIPS_INITIAL(ch_input_for_bips, bips_root_path_obj)

                    def batch_id = "bips_batch_run" // A fixed ID for this batch run
					ch_bips_oligos_csv_result = RUN_BIPS_INITIAL.out.oligos_csv.map { single_oligos_file  -> tuple(batch_id, single_oligos_file) }
					ch_bips_barcoded_csv_result = RUN_BIPS_INITIAL.out.barcoded_csv.map { single_barcode_file  -> tuple(batch_id, single_barcode_file) }
			}

            // ========================
            //   PREPARE INPUT FOR BIPS_CSV_TO_FASTA PROCESS
            // ========================
            // Create a single channel that will feed into BIPS_CSV_TO_FASTA.
            // It will contain data from EITHER bips_then_dolphyn mode OR dolphyn_oligo_only mode.
            if (params.mode == "bips_then_dolphyn" || params.mode == "dolphyn_oligo_only"){
                ch_input_for_csv_conversion = Channel.empty()
                ch_input_for_csv_conversion = ch_input_for_csv_conversion.mix(
                                                ch_bips_oligos_csv_result,      // Has data in bips_then_dolphyn mode
                                                ch_oligos_csv_for_dolphyn_only   // Has data in dolphyn_oligo_only (with CSV input)
                                            )
                // Call the conversion process only ONCE with the merged input channel.
                BIPS_CSV_TO_FASTA(ch_input_for_csv_conversion)

                def ch_fasta_from_csv_conversion = BIPS_CSV_TO_FASTA.out.oligos_fasta

                ch_fasta_for_core_dolphyn = ch_fasta_for_core_dolphyn
                                            .mix(ch_oligos_fasta_for_dolphyn_only)
                                            .mix(ch_fasta_from_csv_conversion)
            }
            

            // ========================
			//   Dolphyn Standalone Execution (on protein FASTA)
			// ========================
			if (params.mode == "dolphyn_standalone") {
					RUN_DOLPHYN_STANDALONE_PREP(ch_input_for_dolphyn_standalone)		// dolphyn_training_data_path_obj // if action_run_dolphyn_standalone needs it)
					ch_dolphyn_standalone_json_result = RUN_DOLPHYN_STANDALONE_PREP.out.dolphyn_json
			}

			if (params.mode == "bips_then_dolphyn" || params.mode == "dolphyn_oligo_only") {
					// If ch_fasta_ready_for_dolphyn is empty, RUN_DOLPHYN_PREDICTION won't run.
					RUN_DOLPHYN_PREDICTION(ch_fasta_for_core_dolphyn)
					ch_dolphyn_json_from_core_prediction  = RUN_DOLPHYN_PREDICTION.out.dolphyn_json

					if (params.mode == "bips_then_dolphyn") {
							// These will only run if RUN_DOLPHYN_PREDICTION ran (i.e., ch_dolphyn_json_result has data)
							DOLPHYN_JSON_TO_CSV(ch_dolphyn_json_from_core_prediction)
							ch_dolphyn_epitope_csv_result = DOLPHYN_JSON_TO_CSV.out.epitope_csv

							// Ensure both channels for join have data
							ch_bips_oligos_csv_result
									.join(ch_dolphyn_epitope_csv_result)
									.set { ch_for_selecting_oligos }
							SELECT_EPITOPE_POSITIVE_BIPS_OLIGOS(ch_for_selecting_oligos)
							ch_selected_bips_oligos_result = SELECT_EPITOPE_POSITIVE_BIPS_OLIGOS.out.selected_bips_oligos_csv

							ch_bips_barcoded_csv_result
									.join(ch_selected_bips_oligos_result)
									.set { ch_for_filtering_barcodes }
							FILTER_BIPS_BARCODES(ch_for_filtering_barcodes)
							ch_final_filtered_barcodes_result = FILTER_BIPS_BARCODES.out.filtered_barcoded_csv
					}
			} else if (params.mode == "bips_then_dolphyn" || params.mode == "dolphyn_oligo_only") {
					log.info "No FASTA input for Dolphyn (channel was empty), skipping Dolphyn steps."
			}

    emit:
        // BIPS general outputs (from bips_only or bips_then_dolphyn)
        bips_oligos_sequence_csv = ch_bips_oligos_csv_result
        bips_barcoded_nuc_csv = ch_bips_barcoded_csv_result
    

        // Dolphyn general output (JSON from oligo prediction)
        dolphyn_oligo_prediction_json = ch_dolphyn_json_from_core_prediction

        // Dolphyn standalone output (JSON from direct protein prediction)
        dolphyn_standalone_prediction_json = ch_dolphyn_standalone_json_result

        // Outputs specific to "bips_then_dolphyn" mode
        intermediate_epitope_list_csv = ch_dolphyn_epitope_csv_result
        selected_bips_oligos_for_barcoding = ch_selected_bips_oligos_result
        final_filtered_barcodes = ch_final_filtered_barcodes_result

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
