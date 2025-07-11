/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { FASTQC                 } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// include { paramsSummaryMap       } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_wookflow_pipeline'
// include { VIRALDB                } from '../subworkflows/local/viraldb/main.nf'
// include { PHIPPERY                } from '../subworkflows/local/phippery/main.nf'
// include { PHIPPERYTOAVARDA        } from '../subworkflows/local/phipperyToAvarda/main.nf'
// include { AVARDA                  } from '../subworkflows/local/AVARDA/main.nf'
// include { POSTAVARDA_WORKFLOW     } from '../subworkflows/local/PostAVARDA/main.nf'

include { BIPS_THEN_DOLPHYN     } from '../subworkflows/local/bipsThenDolphyn/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
log.info"""

#################################
#      Welcome to WookScan      #
#################################

A custom built VirScan analysis suite.

Contributors:
    - Preston (preston@unsw.edu.au)
    - Bea (b.delgado_corrales@unsw.edu.au)
    - Legana (l.fingerhut@unsw.edu.au)
    - Shouyu (Coco) Wei (shouyu.wei@student.unsw.edu.au)

"""


workflow WOOKFLOW {

    // take:
    // ch_samplesheet // channel: samplesheet read in from --input


    log.info "[WORKFLOW] Starting main workflow..."
    log.info "[WORKFLOW] BIPS/Dolphyn mode: ${params.mode ?: 'not set'}"
    log.info "[WORKFLOW] Run Phippery flag: ${params.run_phippery}"
    log.info "[WORKFLOW] Run AVARDA flag: ${params.run_AVARDA}"

    main:
    // ========================================================================================
    //      PRE-PHIP-SEQ: BIPS/DOLPHYN MODULE LOGIC
    // ========================================================================================
    if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only" || params.mode == "dolphyn_standalone" || params.mode == "dolphyn_only") {
        log.info "[WORKFLOW] Entering BIPS/Dolphyn branch based on mode: '${params.mode}'"

        ch_viral_seqs_for_bips = Channel.empty()
        ch_viral_seqs_for_dolphyn_standalone = Channel.empty()
        ch_oligos_fasta_for_dolphyn_only = Channel.empty()
        ch_oligos_csv_for_dolphyn_only_conversion = Channel.empty()


        // --- Prepare Fixed Path Inputs for Processes ---
        def bips_root_path_obj = file(params.bips_root_dir)
        if (!bips_root_path_obj.exists() || !bips_root_path_obj.isDirectory()) {
            error "BIPS root directory specified by params.bips_root_dir not found or not a directory: ${params.bips_root_dir}"
        }

        def helper_script_path_obj = file(params.helper_script)
        if (!helper_script_path_obj.exists()) { error "Helper script not found: ${params.helper_script}" }


        log.info """
            --------------------------------------
            WookScan uses: BuildPhIPSeqLibrary!
            --------------------------------------
            
            BuildPhIPSeqLibrary is developed by:
            kalkairis Iris Kalka, sigallev
            Repository: https://github.com/kalkairis/BuildPhIPSeqLibrary.git            
            ================================
            Input directory  : BuildPhIPSeqLibrary/Input
            Output directory : BuildPhIPSeqLibrary/Output
            
            --------------------------------------
            WookScan Custom Scripts: bipsThenDolphyn
            --------------------------------------
            BIPS output to Dolphyn is developed by:
            Shouyu Wei
            ================================
            mode             : $params.mode
            input directory  : $params.input_viral_seqs_dir
            Output directory : $params.outdir

        """.stripIndent()

        // == Populate Initial Input Channels Based on Mode ==
        if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only") {
            if (params.input_viral_seqs_dir) {
                log.info "Mode '${params.mode}': Attempting Channel.fromPath with params.input_viral_seqs_dir: '${params.input_viral_seqs_dir}'"

                // Create a glob pattern to match all .fa, .fasta, and .csv files in the directory
                def input_glob = "${params.input_viral_seqs_dir}/*.{fa,csv}"
                log.info "Using glob pattern: '${input_glob}'"

                // Channel.fromPath with a glob pattern will emit EACH matching file
                ch_viral_seqs_for_bips = Channel.fromPath(input_glob)
                    // SIMPLIFIED ifEmpty block:
                    .ifEmpty { error "EMPTY CHANNEL: No files (.fa, .fasta, .csv) found in directory provided by --input_viral_seqs_dir: '${params.input_viral_seqs_dir}'. Please check if the file exists and the path is correct relative to Nextflow's launch directory." }
                    .collect() // <-- Collect all emitted files into a single list
                log.info "Channel ch_viral_seqs_for_bips created. It will emit a single list of files."

            } else {
                error "--input_viral_seqs_dir is required for mode '${params.mode}'"
            }
        }


        log.info """
            --------------------------------------
            WookScan uses: Dolphyn!
            --------------------------------------
            
            Dolphyn is developed by:
            Liebhoff, AM. et al
            Repository: https://github.com/kepsi/Dolphyn.git            
            ================================
            Input file  : viral protein sequence file in .fa single line format
            Output file : the predicted epitopes in .json file
            
            --------------------------------------
            WookScan Custom Scripts: bipsThenDolphyn
            --------------------------------------
            The source code of Dolphyn is modified by:
            Shouyu Wei
            
            Modified code includes:
            initEpiPredictor(), getPEDSTrainingSet(), kmer_features_of_protein(), findEpitopes(), saveGlobalEpitopes()
            ================================
            mode             : $params.mode
            input directory  : $params.input_protein_fasta_dir
            Output directory : $params.outdir

        """.stripIndent()

        if (params.mode == "dolphyn_standalone") {
            if (!params.input_protein_fasta_dir) {
                error "Parameter `--input_protein_fasta_dir` (a directory path) must be specified for mode 'dolphyn_standalone'"
            }

            if (params.input_protein_fasta_dir) {
                def input_path = "${params.input_protein_fasta_dir}/*.fa"

                log.info "Mode 'dolphyn_standalone': Creating channel from path/glob: '${input_path}'"
                ch_viral_seqs_for_dolphyn_standalone = Channel.fromPath(input_path)
                    // SIMPLIFIED ifEmpty block:
                    .ifEmpty { error "EMPTY CHANNEL: Channel.fromPath (for dolphyn_standalone) created an empty channel. Path/Pattern was: '${input_path}'. Please check existence and path." }
                    .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_viral_seqs_for_dolphyn_standalone created."
            }
        }


        log.info """
            --------------------------------------
            WookScan uses: Dolphyn!
            --------------------------------------
            
            Dolphyn is developed by:
            Liebhoff, AM. et al
            Repository: https://github.com/kepsi/Dolphyn.git            
            ================================
            Input file  : oligo sequence file in .fa or .csv single line format
            Output file : the predicted epitopes in .json file
            
            --------------------------------------
            WookScan Custom Scripts: bipsThenDolphyn
            --------------------------------------
            The source code of Dolphyn is modified by:
            Shouyu Wei
            
            Modified code includes:
            initEpiPredictor(), getPEDSTrainingSet(), kmer_features_of_protein(), findEpitopes(), saveGlobalEpitopes()
            ================================
            mode             : $params.mode
            input directory  : $params.input_oligos_fasta_dir OR $params.input_oligos_csv_dir
            Output directory : $params.outdir

        """.stripIndent()

        if (params.mode == "dolphyn_only") {
            // Check that only one of the two possible directory inputs is provided
            if (params.input_oligos_fasta_dir && params.input_oligos_csv_dir) {
                error "For 'dolphyn_only' mode, please specify only ONE of --input_oligos_fasta_dir or --input_oligos_csv_dir."
            }
            if (!params.input_oligos_fasta_dir && !params.input_oligos_csv_dir) {
                error "For 'dolphyn_only' mode, either --input_oligos_fasta_dir or --input_oligos_csv_dir must be specified."
            }
            


            if (params.input_oligos_fasta_dir) {
                def input_path = "${params.input_oligos_fasta_dir}/*.fa"
                log.info "Mode 'dolphyn_only': Creating channel from oligo FASTA path/glob: '${input_path}'"
                
                ch_oligos_fasta_for_dolphyn_only = Channel.fromPath(input_path)
                                                    .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for oligos_fasta) created an empty channel. Path/Pattern was: '${input_path}'." }
                                                    .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_oligos_fasta_for_dolphyn_only  created."

            } else if (params.input_oligos_csv_dir) {
                def input_path = params.input_bips_oligos_csv ?: "${params.input_oligos_csv_dir}/*.csv"
                log.info "Mode 'dolphyn_only': Creating channel from BIPS CSV path/glob: '${input_path}'"
                ch_oligos_csv_for_dolphyn_only_conversion = Channel.fromPath(input_path)
                                                        .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for bips_csv) created an empty channel. Path/Pattern was: '${input_path}'." }
                                                        .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_oligos_csv_for_dolphyn_only_conversion created."
            }
        }

        // --- Call bipsThenDolphyn subworkflow ---
        // Only call it if the mode matches one of the ones it handles
        if (params.mode.startsWith("bips") || params.mode.startsWith("dolphyn")) {
            log.info "Starting BIPS/Dolphyn module..."
            BIPS_THEN_DOLPHYN (
                ch_viral_seqs_for_bips,
                ch_viral_seqs_for_dolphyn_standalone,
                ch_oligos_fasta_for_dolphyn_only,
                ch_oligos_csv_for_dolphyn_only_conversion,
                bips_root_path_obj,
                helper_script_path_obj
            )

            // You can use the outputs from your subworkflow here if needed
            BIPS_THEN_DOLPHYN.out.bips_oligos_sequence_csv.view { "Generated BIPS oligos CSV: $it" }
        }
    }
    



    // ========================================================================================
    //      POST-PHIP-SEQ: PHIPPERY/AVARDA MODULE LOGIC
    // ========================================================================================

    // ch_versions = Channel.empty()
    // ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    
    // FASTQC (
    //     ch_samplesheet
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
    //     .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    
    // ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    // ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    // summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )

    // emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    // versions       = ch_versions                 // channel: [ path(versions.yml) ]


    //
    // RUN: ViralDB - Build Viral Library support files for AVARDA
    //
    
    // TODO


    //
    // RUN: fastp to clean files
    //
    
    // TODO:
    // Need to figure out how to connect the fastq files to fastp. 
    // After fastp, we need to manipulat the sample table information so that
    // phippery can read the trimmed fastq output from fastp.


    //
    // RUN: PHIPPERY
    //    
//     if(params.run_phippery){        
//         log.info """
//             --------------------------------------
//             WookScan uses: P H I P - F L O W!
//             --------------------------------------
            
//             Phippery & phip-flow is developed by:
//             Matsen, Overbaugh, and Minot Labs
//             Fred Hutchinson CRC, Seattle WA
//             Repository: https://github.com/matsengrp/phip-flow            
//             ================================
//             sample_table    : $params.sample_table
//             peptide_table   : $params.peptide_table
//             results         : $params.results
//             reads_prefix    : $params.reads_prefix
//             publishDir      : $params.results
            
//             --------------------------------------
//             WookScan Custom Scripts: PhipOut
//             --------------------------------------
//             Phippery output to AVARDA is developed by:
//             Preston Leung
//             ================================
//             user_pep_id     : $params.user_pep_id
//             publishDir      : $params.results

//         """.stripIndent()
//         PHIPPERY()        
//         PHIPPERYTOAVARDA(PHIPPERY.out)
//     }

//     if(params.run_AVARDA){        
//         if(params.run_phippery){
//             // If we take stuff directly from phippery output
//             log.info """
//             --------------------------------------
//             WookScan uses: AVARDA: PHIPPERY-AVARDA
//             --------------------------------------
//             AVARDA is developed by:
//             Monaco et al.
//             Repository: https://github.com/drmonaco/AVARDA

//             Modification performed by:
//             Preston Leung
//             ================================
//             virlib      : From PHIPPERYTOAVARDA.out.virlib
//             edgeRhits   : From PHIPPERYTOAVARDA.out.edgeRhits
//             publishDir  : $params.out_path
            
//             """.stripIndent()
//             virlib = PHIPPERYTOAVARDA.out.virlib
//             edgeRhits = PHIPPERYTOAVARDA.out.edgeRhits
//         }else{
//             // We're running AVARDA by itself
//             log.info """
//             --------------------------------------
//             WookScan uses: AVARDA: AVARDA ONLY
//             --------------------------------------
//             AVARDA is developed by:
//             Monaco et al.

//             Modification performed by:
//             Preston Leung
//             ================================
//             virlib      : $params.avarda_names
//             edgeRhits   : $params.case_path
//             publishDir  : $params.out_path

//             """.stripIndent()
//             virlib = Channel.fromPath(params.avarda_names)
//             edgeRhits =  Channel.fromPath(params.case_path)
//         }        
//         AVARDA(virlib, edgeRhits)
//         POSTAVARDA_WORKFLOW()
//     }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
