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
*   Dolphyn-related: dolphyn_standalone, dolphyn_only.
*   
*/

nextflow.enable.dsl=2


// --- Workflow ---
workflow BIPS_THEN_DOLPHYN {
    take:
        ch_input_for_bips  // ch_viral_seqs_for_bips		// For modes: bips_then_dolphyn, bips_only
        ch_input_for_dolphyn_standalone // ch_viral_seqs_for_dolphyn_standalone		// For mode: dolphyn_standalone
        ch_oligos_fasta_for_dolphyn_only		// For mode: dolphyn_only (oligo FASTA)
        ch_bips_csv_for_dolphyn_only		// For mode: dolphyn_only (oligo CSV)


        //ch_input_files
        bips_root_path_obj // Path object for BIPS vendored code
        helper_script_path_obj // Path object for your Python helper



    main:
    log.info("${params.mode}")
    log.info("${ch_oligos_fasta_for_dolphyn_only}")
    log.info("${ch_bips_csv_for_dolphyn_only}")
    // == Initialize Core Channels ==
    // These will be populated based on the pipeline mode.
    // ch_viral_seqs_for_bips = Channel.empty()
    // ch_viral_seqs_for_dolphyn_standalone = Channel.empty()
    // ch_oligos_fasta_for_dolphyn_only = Channel.empty()
    // ch_bips_csv_for_dolphyn_only_conversion = Channel.empty()



    // --- Prepare Fixed Path Inputs for Processes ---
    // def bips_root_path_obj = file(params.bips_root_dir)
    // if (!bips_root_path_obj.exists() || !bips_root_path_obj.isDirectory()) {
    //     error "BIPS root directory specified by params.bips_root_dir not found or not a directory: ${params.bips_root_dir}"
    // }

    // def helper_script_path_obj = file(params.helper_script)
    // if (!helper_script_path_obj.exists()) { error "Helper script not found: ${params.helper_script}" }

    // // == Populate Initial Input Channels Based on Mode ==
    // if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only") {
    //     if (params.input_viral_seqs) {
    //         log.info "Mode '${params.mode}': Attempting Channel.fromPath with params.input_viral_seqs: '${params.input_viral_seqs}'"
    //         ch_viral_seqs_for_bips = Channel.fromPath(params.input_viral_seqs)
    //             // SIMPLIFIED ifEmpty block:
    //             .ifEmpty { error "EMPTY CHANNEL: Channel.fromPath (for viral_seqs in mode '${params.mode}') created an empty channel. Input parameter was: '${params.input_viral_seqs}'. Please check if the file exists and the path is correct relative to Nextflow's launch directory." }
    //             .map { f -> tuple(f.baseName, f) }
    //         log.info "Channel ch_viral_seqs_for_bips potentially created."

    //     } else {
    //         error "--input_viral_seqs is required for mode '${params.mode}'"
    //     }
    // }

    // if (params.mode == "dolphyn_standalone") {
    //     if (params.input_viral_seqs) {
    //         log.info "Mode 'dolphyn_standalone': Attempting Channel.fromPath with params.input_viral_seqs: '${params.input_viral_seqs}'"
    //         ch_viral_seqs_for_dolphyn_standalone = Channel.fromPath(params.input_viral_seqs)
    //             // SIMPLIFIED ifEmpty block:
    //             .ifEmpty { error "EMPTY CHANNEL: Channel.fromPath (for dolphyn_standalone) created an empty channel. Input parameter was: '${params.input_viral_seqs}'. Please check file existence and path." }
    //             .map { f -> tuple(f.baseName, f) }
    //         log.info "Channel ch_viral_seqs_for_dolphyn_standalone potentially created."
    //     }
    // }

    // if (params.mode == "dolphyn_only") {
    //     if (params.input_oligos_fasta) {
    //         log.info "Mode 'dolphyn_only': Attempting Channel.fromPath with params.input_oligos_fasta: '${params.input_oligos_fasta}'"
    //         ch_oligos_fasta_for_dolphyn_only = Channel.fromPath(params.input_oligos_fasta)
    //                                               .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for oligos_fasta) created an empty channel. Input parameter was: '${params.input_oligos_fasta}'." }
    //                                               .map { f -> tuple(f.baseName, f) }
    //         log.info "Channel ch_oligos_fasta_for_dolphyn_only potentially created."

    //     } else if (params.input_bips_oligos_csv) {
    //         log.info "Mode 'dolphyn_only': Attempting Channel.fromPath with params.input_bips_oligos_csv: '${params.input_bips_oligos_csv}'"
    //         ch_bips_csv_for_dolphyn_only_conversion = Channel.fromPath(params.input_bips_oligos_csv)
    //                                                       .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for bips_csv) created an empty channel. Input parameter was: '${params.input_bips_oligos_csv}'." }
    //                                                       .map { f -> tuple(f.baseName, f) }
    //         log.info "Channel ch_bips_csv_for_dolphyn_only_conversion potentially created."
    //     }
    // }



			// == Initialize Core Channels ==
			// These will be populated based on the pipeline mode.
			ch_bips_oligos_csv_result = Channel.empty()
			ch_bips_barcoded_csv_result = Channel.empty()

            ch_bips_full_output_dir = Channel.empty()

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
                    def batch_id = "bips_batch_run" // A fixed ID for this batch run

					RUN_BIPS_INITIAL(ch_input_for_bips, bips_root_path_obj)
					ch_bips_oligos_csv_result = RUN_BIPS_INITIAL.out.oligos_csv.map { file -> tuple(batch_id, file) }
					ch_bips_barcoded_csv_result = RUN_BIPS_INITIAL.out.barcoded_csv.map { file -> tuple(batch_id, file) }
                    ch_bips_full_output_dir = RUN_BIPS_INITIAL.out.bips_full_outputs_dir.map { file -> tuple(batch_id, file) }

                    PUBLISH_BIPS_RESULTS(ch_bips_full_output_dir)
                    PUBLISH_BIPS_RESULTS.out.view { "BIPS publishing signal file: $it" }

					// if (params.mode == "bips_then_dolphyn") {
					// 		BIPS_CSV_TO_FASTA(ch_bips_oligos_csv_result)
					// 		ch_fasta_for_core_dolphyn = BIPS_CSV_TO_FASTA.out.oligos_fasta
					// }
			}

            // ========================
            //   PREPARE INPUT FOR BIPS_CSV_TO_FASTA PROCESS
            // ========================
            // Create a single channel that will feed into BIPS_CSV_TO_FASTA.
            // It will contain data from EITHER bips_then_dolphyn mode OR dolphyn_only mode.
            ch_input_for_csv_conversion = Channel.empty()
            ch_input_for_csv_conversion = ch_input_for_csv_conversion.mix(
                                            ch_bips_oligos_csv_result,      // Has data in bips_then_dolphyn mode
                                            ch_bips_csv_for_dolphyn_only   // Has data in dolphyn_only (with CSV input)
                                        )
            // Call the conversion process only ONCE with the merged input channel.
            BIPS_CSV_TO_FASTA(ch_input_for_csv_conversion)
			

			// ========================
			//   Dolphyn Direct Input Prep
			// ========================
            // We don't need an `if (params.mode == 'dolphyn_only')` block here anymore,
            // because the `when:` directive on the process handles the logic.

            // This process call will ONLY be activated if the 'when:' condition inside it is met.
            // We pass the channel containing the user-provided BIPS CSV. If the channel is empty
            // (because the user didn't provide the param), this process won't run.
            //BIPS_CSV_TO_FASTA(ch_bips_csv_for_dolphyn_only)
            def ch_fasta_from_csv_conversion = BIPS_CSV_TO_FASTA.out.oligos_fasta

            ch_fasta_for_core_dolphyn = ch_fasta_for_core_dolphyn
                                        .mix(ch_oligos_fasta_for_dolphyn_only)
                                        .mix(ch_fasta_from_csv_conversion)




            // ========================
			//   Dolphyn Standalone Execution (on protein FASTA)
			// ========================
			if (params.mode == "dolphyn_standalone") {
					RUN_DOLPHYN_STANDALONE_PREP(ch_input_for_dolphyn_standalone)		// dolphyn_training_data_path_obj // if action_run_dolphyn_standalone needs it)
					ch_dolphyn_standalone_json_result = RUN_DOLPHYN_STANDALONE_PREP.out.dolphyn_json
			}



			// if (params.mode == "dolphyn_only") {
            //     //def ch_bips_csv_for_dolphyn_only = 

            //     //use 'when' here to seperate ch_oligos_fasta_for_dolphyn_only and ch_bips_csv_for_dolphyn_only!
            //     //channel.path = 'assets/dol_test_data/oligos_sequence.csv'

            //     if (ch_oligos_fasta_for_dolphyn_only.ifEmpty(false)) {   // Ask
            //             ch_fasta_for_core_dolphyn = ch_oligos_fasta_for_dolphyn_only
            //     }

            //     if (ch_bips_csv_for_dolphyn_only.ifEmpty(false)) {
            //         // Convert user-provided BIPS CSV to FASTA for Dolphyn
            //         BIPS_CSV_TO_FASTA(ch_bips_csv_for_dolphyn_only)
            //         ch_fasta_for_core_dolphyn = BIPS_CSV_TO_FASTA.out.oligos_fasta
            //     } else {
            //         log.info "Error in dolphyn_only."
            //     }
                    
			// }

            // Directly connect the BIPS CSV to FASTA conversion process
            // BIPS_CSV_TO_FASTA(ch_bips_csv_for_dolphyn_only)
            // def ch_fasta_from_bips_csv = BIPS_CSV_TO_FASTA.out.oligos_fasta

            // // Mix all potential sources for the core Dolphyn input.
            // // Only the channels that actually receive data from the main workflow will contribute.
            // ch_fasta_for_core_dolphyn = ch_fasta_for_core_dolphyn
            //                             .mix(BIPS_CSV_TO_FASTA.out.oligos_fasta) // Add results from bips_then_dolphyn
            //                             .mix(ch_oligos_fasta_for_dolphyn_only) // Add results from dolphyn_only (direct fasta)
            //                             .mix(ch_fasta_from_bips_csv) // Add results from dolphyn_only (csv conversion)

			// ========================
			//    Dolphyn Execution & Downstream (if applicable)
			// ========================
			// Check if ch_fasta_ready_for_dolphyn actually received data
			// Use .take(1) to ensure it proceeds if there's at least one item,
			// then filter it back if needed, or use a flag.
			// A simpler check for "will this branch run":
			// ch_run_dolphyn_signal = ch_fasta_ready_for_dolphyn
			//     .count() // Counts items, emits the count
			//     .map { count -> count > 0 }
			//boolean run_dolphyn_branch = (params.mode == "bips_then_dolphyn" || params.mode == "dolphyn_only") && !ch_fasta_ready_for_dolphyn.isEmpty().toBlocking().first()

			if (params.mode == "bips_then_dolphyn" || params.mode == "dolphyn_only") {
					ch_fasta_for_core_dolphyn
							.ifEmpty { log.info "No oligo FASTA input for core Dolphyn, skipping RUN_DOLPHYN_PREDICTION and subsequent steps." }
							.set { ch_valid_oligo_fasta_for_dolphyn }


					// If ch_fasta_ready_for_dolphyn is empty, RUN_DOLPHYN_PREDICTION won't run.
					RUN_DOLPHYN_PREDICTION(ch_valid_oligo_fasta_for_dolphyn)
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
			} else if (params.mode == "bips_then_dolphyn" || params.mode == "dolphyn_only") {
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

// --- Processes ---

process RUN_BIPS_INITIAL {
    // The tag now represents a batch run, so sample_id is not directly applicable.
    // We can use a fixed tag or one based on the mode.
    tag "BIPS Viral batch run for mode ${params.mode}"

    input:
    //tuple val(sample_id), path(viral_seq_file)
    path file_list
    path bips_code_dir

    output:
    // These paths are now relative to the Nextflow work directory after being copied from BIPS_run
    // They are a single set of results for the entire batch.
    path(params.bips_oligos_sequence_csv_name), emit: oligos_csv
    path(params.bips_barcoded_nuc_file_csv_name), emit: barcoded_csv
    path ("${bips_code_dir.baseName}/Data/Output"), emit: bips_full_outputs_dir

    script:
    def bips_staged_name = bips_code_dir.baseName // e.g., "BuildPhIPSeqLibrary"
    def oligos_out_name = params.bips_oligos_sequence_csv_name
    def barcoded_out_name = params.bips_barcoded_nuc_file_csv_name

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

    # Check if the BIPS output directory itself was created
    if [ ! -d "${bips_staged_name}/Data/Output" ]; then
        echo "ERROR: BIPS did not create the output directory: ${bips_staged_name}/Data/Output" >&2
        exit 1
    fi

    if [ -f "${bips_staged_name}/Data/Output/${oligos_out_name}" ]; then
        cp "${bips_staged_name}/Data/Output/${oligos_out_name}" .  # Copy to current dir (root)
        echo "Copied ${oligos_out_name} to work dir root for output."
    else
        echo "ERROR: BIPS output ${oligos_out_name} not found!" >&2
        exit 1
    fi

    if [ -f "${bips_staged_name}/Data/Output/${barcoded_out_name}" ]; then
        cp "${bips_staged_name}/Data/Output/${barcoded_out_name}" .  # Copy to current dir (root)
        echo "Copied ${barcoded_out_name} to work dir root for output."
    else
        echo "ERROR: BIPS output ${barcoded_out_name} not found!" >&2
        exit 1
    fi
    """
}

process PUBLISH_BIPS_RESULTS {
    // This is a simple utility process, tag can be descriptive
    tag "Publishing BIPS full results for $sample_id"

    // No publishDir needed here, as its whole job is to publish!
    // This process runs locally and doesn't need a complex environment.
    // It just needs 'cp'.

    input:
    // It receives the sample_id and the Path object for the directory to be published
    tuple val(sample_id), path(bips_output_dir)

    output:
    // This process doesn't need to produce an output for other processes,
    // but we can emit a signal file to indicate completion.
    path "${sample_id}.bips_published.txt"

    script:
    // Define the final target directory path
    def final_path = "${params.outdir}/${sample_id}/bips_run_outputs"
    """
    echo "Publishing BIPS full results for sample ${sample_id}"
    echo "Source directory (from channel): ${bips_output_dir}"
    echo "Target directory (final output): ${final_path}"

    # Create the target directory
    mkdir -p "${final_path}"

    # Copy the *contents* of the source directory to the target directory
    # The '/*' at the end is important for copying contents
    cp -rL "${bips_output_dir}"/* "${final_path}/"

    # Create a small file to signal that publishing is done
    echo "Published on \$(date)" > "${sample_id}.bips_published.txt"
    """
}

process BIPS_CSV_TO_FASTA {
    tag "$sample_id"
    publishDir "${params.outdir}/dolphyn_prep", mode: 'copy', pattern: "*.fasta"

    when:
        (params.mode == 'bips_then_dolphyn') || (params.mode == 'dolphyn_only' && params.input_bips_oligos_csv)
        // && params.bips_oligos_sequence_csv_name

    input:
    tuple val(sample_id), path(bips_oligos_csv) // from RUN_BIPS_INITIAL

    output:
    tuple val(sample_id), path("${sample_id}.oligos_for_dolphyn.fasta"), emit: oligos_fasta

    script:
    """
    ${params.helper_script} bips_csv_to_fasta \\
        --bips_oligos_csv ${bips_oligos_csv} \\
        --output_fasta ${sample_id}.oligos_for_dolphyn.fasta
    """
}

process RUN_DOLPHYN_PREDICTION {
    conda "bips_environment.yml"
    
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
    conda "bips_environment.yml"


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
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
