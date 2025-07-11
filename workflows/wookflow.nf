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
include { PHIPPERY                } from '../subworkflows/local/phippery/main.nf'
include { PHIPPERYTOAVARDA        } from '../subworkflows/local/phipperyToAvarda/main.nf'
include { AVARDA                  } from '../subworkflows/local/AVARDA/main.nf'
include { POSTAVARDA_WORKFLOW     } from '../subworkflows/local/PostAVARDA/main.nf'

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

    main:
    // ========================================================================================
    //      PRE-PHIP-SEQ: BIPS/DOLPHYN MODULE LOGIC
    // ========================================================================================
        ch_viral_seqs_for_bips = Channel.empty()
        ch_viral_seqs_for_dolphyn_standalone = Channel.empty()
        ch_oligos_fasta_for_dolphyn_only = Channel.empty()
        ch_bips_csv_for_dolphyn_only_conversion = Channel.empty()



        // --- Prepare Fixed Path Inputs for Processes ---
        def bips_root_path_obj = file(params.bips_root_dir)
        if (!bips_root_path_obj.exists() || !bips_root_path_obj.isDirectory()) {
            error "BIPS root directory specified by params.bips_root_dir not found or not a directory: ${params.bips_root_dir}"
        }

        def helper_script_path_obj = file(params.helper_script)
        if (!helper_script_path_obj.exists()) { error "Helper script not found: ${params.helper_script}" }

        // == Populate Initial Input Channels Based on Mode ==
        if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only") {
            if (params.input_viral_seqs) {
                log.info "Mode '${params.mode}': Attempting Channel.fromPath with params.input_viral_seqs: '${params.input_viral_seqs}'"
                ch_viral_seqs_for_bips = Channel.fromPath(params.input_viral_seqs)
                    // SIMPLIFIED ifEmpty block:
                    .ifEmpty { error "EMPTY CHANNEL: Channel.fromPath (for viral_seqs in mode '${params.mode}') created an empty channel. Input parameter was: '${params.input_viral_seqs}'. Please check if the file exists and the path is correct relative to Nextflow's launch directory." }
                    .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_viral_seqs_for_bips potentially created."

            } else {
                error "--input_viral_seqs is required for mode '${params.mode}'"
            }
        }

        if (params.mode == "dolphyn_standalone") {
            if (params.input_viral_seqs) {
                log.info "Mode 'dolphyn_standalone': Attempting Channel.fromPath with params.input_viral_seqs: '${params.input_viral_seqs}'"
                ch_viral_seqs_for_dolphyn_standalone = Channel.fromPath(params.input_viral_seqs)
                    // SIMPLIFIED ifEmpty block:
                    .ifEmpty { error "EMPTY CHANNEL: Channel.fromPath (for dolphyn_standalone) created an empty channel. Input parameter was: '${params.input_viral_seqs}'. Please check file existence and path." }
                    .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_viral_seqs_for_dolphyn_standalone potentially created."
            }
        }

        if (params.mode == "dolphyn_only") {
            if (params.input_oligos_fasta) {
                log.info "Mode 'dolphyn_only': Attempting Channel.fromPath with params.input_oligos_fasta: '${params.input_oligos_fasta}'"
                ch_oligos_fasta_for_dolphyn_only = Channel.fromPath(params.input_oligos_fasta)
                                                    .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for oligos_fasta) created an empty channel. Input parameter was: '${params.input_oligos_fasta}'." }
                                                    .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_oligos_fasta_for_dolphyn_only potentially created."

            } else if (params.input_bips_oligos_csv) {
                log.info "Mode 'dolphyn_only': Attempting Channel.fromPath with params.input_bips_oligos_csv: '${params.input_bips_oligos_csv}'"
                ch_bips_csv_for_dolphyn_only_conversion = Channel.fromPath(params.input_bips_oligos_csv)
                                                            .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for bips_csv) created an empty channel. Input parameter was: '${params.input_bips_oligos_csv}'." }
                                                            .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_bips_csv_for_dolphyn_only_conversion potentially created."
            }
        }

        // --- Call bipsThenDolphyn subworkflow ---
        // Only call it if the mode matches one of the ones it handles
        if (params.mode.startsWith("bips") || params.mode.startsWith("dolphyn")) {
            log.info "Starting BIPS/Dolphyn module..."
            BIPS_THEN_DOLPHYN_SUB (
                ch_viral_seqs_for_bips,
                ch_viral_seqs_for_dolphyn_standalone,
                ch_oligos_fasta_for_dolphyn_only,
                ch_bips_csv_for_dolphyn_only_conversion,
                bips_root_path_obj,
                helper_script_path_obj
            )

            // You can use the outputs from your subworkflow here if needed
            BIPS_THEN_DOLPHYN.out.bips_oligos_sequence_csv.view { "Generated BIPS oligos CSV: $it" }
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
    if(params.run_phippery){        
        log.info """
            --------------------------------------
            WookScan uses: P H I P - F L O W!
            --------------------------------------
            
            Phippery & phip-flow is developed by:
            Matsen, Overbaugh, and Minot Labs
            Fred Hutchinson CRC, Seattle WA
            Repository: https://github.com/matsengrp/phip-flow            
            ================================
            sample_table    : $params.sample_table
            peptide_table   : $params.peptide_table
            results         : $params.results
            reads_prefix    : $params.reads_prefix
            publishDir      : $params.results
            
            --------------------------------------
            WookScan Custom Scripts: PhipOut
            --------------------------------------
            Phippery output to AVARDA is developed by:
            Preston Leung
            ================================
            user_pep_id     : $params.user_pep_id
            publishDir      : $params.results

        """.stripIndent()
        PHIPPERY()        
        PHIPPERYTOAVARDA(PHIPPERY.out)
    }

    if(params.run_AVARDA){        
        if(params.run_phippery){
            // If we take stuff directly from phippery output
            log.info """
            --------------------------------------
            WookScan uses: AVARDA: PHIPPERY-AVARDA
            --------------------------------------
            AVARDA is developed by:
            Monaco et al.
            Repository: https://github.com/drmonaco/AVARDA

            Modification performed by:
            Preston Leung
            ================================
            virlib      : From PHIPPERYTOAVARDA.out.virlib
            edgeRhits   : From PHIPPERYTOAVARDA.out.edgeRhits
            publishDir  : $params.out_path
            
            """.stripIndent()
            virlib = PHIPPERYTOAVARDA.out.virlib
            edgeRhits = PHIPPERYTOAVARDA.out.edgeRhits
        }else{
            // We're running AVARDA by itself
            log.info """
            --------------------------------------
            WookScan uses: AVARDA: AVARDA ONLY
            --------------------------------------
            AVARDA is developed by:
            Monaco et al.

            Modification performed by:
            Preston Leung
            ================================
            virlib      : $params.avarda_names
            edgeRhits   : $params.case_path
            publishDir  : $params.out_path

            """.stripIndent()
            virlib = Channel.fromPath(params.avarda_names)
            edgeRhits =  Channel.fromPath(params.case_path)
        }        
        AVARDA(virlib, edgeRhits)
        POSTAVARDA_WORKFLOW()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
