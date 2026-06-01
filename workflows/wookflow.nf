/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { VIRALDB                } from '../subworkflows/local/viraldb/main.nf'
include { CUTADAPT_WORKFLOW   } from '../subworkflows/local/cutadapt/main.nf'
include { FASTP_WORKFLOW      } from '../subworkflows/local/fastp/main.nf'
include { PEARMERGE_WORKFLOW  } from '../subworkflows/local/pearmerge/main.nf'
include { PHIPPERY            } from '../subworkflows/local/phippery/main.nf'
include { PHIPPERYTOAVARDA    } from '../subworkflows/local/phipperyToAvarda/main.nf'
include { AVARDA              } from '../subworkflows/local/AVARDA/main.nf'
// may remove this later
// include { POSTAVARDA_WORKFLOW     } from '../subworkflows/local/PostAVARDA/main.nf'
include { BIPS_THEN_DOLPHYN } from '../subworkflows/local/bipsThenDolphyn/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// HOW TO USE PROFILES TO DO CRAP 
// FIX THIS NF FILE TO THIS FORMAT

// include { SUB_A } from './subworkflows/analysis_a'
// include { SUB_B } from './subworkflows/analysis_b'

// workflow {
//     // 1. Get a list of active profiles
//     def active_profiles = workflow.profile.tokenize(',')

//     // 2. Conditional logic to run subworkflows
//     if ( active_profiles.contains('analysis_a') ) {
//         SUB_A(ch_input)
//     } 
    
//     if ( active_profiles.contains('analysis_b') ) {
//         SUB_B(ch_input)
//     }
// }



workflow WOOKFLOW {

    // take:
    // ch_samplesheet 
    // channel: samplesheet read in from --input    

    main:
    
    log.info"""

        #################################
        #      Welcome to WookScan      #
        #################################

        A custom built PhIPSeq analysis suite. Currently supports:
            - VirScan

        Contributors:
            - Preston (preston@unsw.edu.au)
            - Bea (b.delgado_corrales@unsw.edu.au)
            - Legana (l.fingerhut@unsw.edu.au)
            - Shouyu (Coco) Wei (shouyu.wei@student.unsw.edu.au)

    """



    def active_profiles = workflow.profile.tokenize(',')
    def defined_profiles = [
        'virscan', 'huscan', 'avarda_only',
        'standard', 'debug','conda',
        'docker','singularity','arm',
        'test','test_full'
    ]
    def invalid_profiles = active_profiles - defined_profiles
    if(invalid_profiles.size() > 0) {
        error "CRITICAL ERROR: Invalid profile(s) specified: ${invalid_profiles}. Please check your command line parameters and ensure all profiles are valid."
    }
    // Some profile checks to avoid conflicting / incompatible profiles
    if (active_profiles.contains('virscan') && active_profiles.contains('avarda_only') ) {
        error "CRITICAL ERROR: Profile 'virscan' includes 'avarda_only' as an option. Please add '--run_AVARDA true' in command line params and remove 'avarda_only' from the profile param."
    }
    if (active_profiles.contains('huscan') && active_profiles.contains('avarda_only')) {
        error "CRITICAL ERROR: Profile 'huscan' is not compatible with avarda tool. Please omit 'avarda' when specifying profiles in command line params."
    }    
    if ( active_profiles.contains('virscan') && active_profiles.contains('huscan') ) {
        error "CRITICAL ERROR: Profiles 'virscan' and 'huscan' are mutually exclusive. Please choose only one."
    }  
    

    log.info "[WORKFLOW] Profile checks completed: Starting main workflow!"
    log.info "[WORKFLOW] Profiles: ${active_profiles}"
    // log.info "[WORKFLOW] BIPS/Dolphyn mode: ${params.mode ?: 'not set'}"
    // log.info "[WORKFLOW] Run Phippery flag: ${params.run_phippery}"
    // log.info "[WORKFLOW] Run AVARDA flag: ${params.run_AVARDA}"
    
    // ==============================================================================================
    //      PRE-PHIP-SEQ: BIPS/DOLPHYN MODULE LOGIC - NEED TO REWORK THIS LATER TO FIT INTO PROFILES
    // ==============================================================================================
    if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only" || params.mode == "dolphyn_standalone" || params.mode == "dolphyn_oligo_only") {
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

        // == Populate Initial Input Channels Based on Mode ==
        if (params.mode == "bips_then_dolphyn" || params.mode == "bips_only") {
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

            if (params.input_viral_seqs_dir) {
                log.info "Mode '${params.mode}': Creating BATCH channel from directory: '${params.input_viral_seqs_dir}'"

                def input_glob = "${params.input_viral_seqs_dir}/*.{fa,csv}"
                log.info "Using path glob pattern: '${input_glob}'"

                // --- THIS IS THE SINGLE, CORRECT LOGIC FOR BATCH MODE ---
                // 1. Create a channel from the glob pattern (emits multiple files).
                // 2. Collect all files into a single list.
                // 3. This channel will now emit ONE item: the list of files.
                ch_viral_seqs_for_bips = Channel.fromPath(input_glob)
                    .ifEmpty { error "EMPTY CHANNEL: No files (.fa, .csv) found in directory: '${params.input_viral_seqs_dir}'." }
                    .collect()

                log.info "Channel ch_viral_seqs_for_bips created. It will emit a single list of files for batch processing."

            } else {
                error "--input_viral_seqs_dir is required for mode '${params.mode}'"
            }
        }

        if (params.mode == "dolphyn_standalone") {
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

        if (params.mode == "dolphyn_oligo_only") {
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
            // Check that only one of the two possible directory inputs is provided
            if (params.input_oligos_fasta_dir && params.input_oligos_csv_dir) {
                error "For 'dolphyn_oligo_only' mode, please specify only ONE of --input_oligos_fasta_dir or --input_oligos_csv_dir."
            }
            if (!params.input_oligos_fasta_dir && !params.input_oligos_csv_dir) {
                error "For 'dolphyn_oligo_only' mode, either --input_oligos_fasta_dir or --input_oligos_csv_dir must be specified."
            }
            if (params.input_oligos_fasta_dir) {
                def input_path = "${params.input_oligos_fasta_dir}/*.fa"
                log.info "Mode 'dolphyn_oligo_only': Creating channel from oligo FASTA path/glob: '${input_path}'"
                
                ch_oligos_fasta_for_dolphyn_only = Channel.fromPath(input_path)
                                                    .ifEmpty{ error "EMPTY CHANNEL: Channel.fromPath (for oligos_fasta) created an empty channel. Path/Pattern was: '${input_path}'." }
                                                    .map { f -> tuple(f.baseName, f) }
                log.info "Channel ch_oligos_fasta_for_dolphyn_only  created."
            } else if (params.input_oligos_csv_dir) {
                // def input_path = params.input_bips_oligos_csv ?: "${params.input_oligos_csv_dir}/*.csv"
                def input_path = "${params.input_oligos_csv_dir}/*.csv"
                log.info "Mode 'dolphyn_oligo_only': Creating channel from BIPS CSV path/glob: '${input_path}'"
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
        }
    }   

    // ========================================================================================
    //      POST-PHIP-SEQ: PHIPPERY/AVARDA MODULE LOGIC
    // ======================================================================================== 

    // VirScan mode: run Phippery and AVARDA (optional) with custom WookScan modifications
    if(params.runtype == 'virscan'){        
        // Print out the message for running phippery            
        log.info """
            --------------------------------------
            WookScan uses: P H I P - F L O W!
            Runtype         : $params.runtype
            --------------------------------------            
            Phippery & phip-flow is developed by:
            -   Matsen, Overbaugh, and Minot Labs
            -   Fred Hutchinson CRC, Seattle WA
            Repository: 
            -   https://github.com/matsengrp/phip-flow            
            ================================
            sample_table    : $params.sample_table
            peptide_table   : $params.peptide_table            
            publishDir      : $params.results
            
            --------------------------------------
            WookScan Custom Scripts: PhipOut
            --------------------------------------
            Phippery output to AVARDA is developed by:
            -   Preston Leung
            ================================
            user_pep_id     : $params.user_pep_id
            publishDir      : $params.results
        """.stripIndent()        
        
        if(params.run_AVARDA){ 
            // If we take stuff directly from phippery output
            log.info """
                --------------------------------------
                WookScan uses   : AVARDA: PHIPPERY-AVARDA
                Runtype         : $params.runtype
                --------------------------------------
                AVARDA is developed by:
                -   Monaco et al.
                Repository: 
                -   https://github.com/drmonaco/AVARDA
                Modification performed by:
                -   Preston Leung
                ================================
                virlib          : From PHIPPERYTOAVARDA.out.virlib
                edgeRhits       : From PHIPPERYTOAVARDA.out.edgeRhits
                publishDir      : $params.out_path            
            """.stripIndent()
        }
        sample_ch = Channel.fromPath(params.sample_table)
        // Running Fastp -> Phippery
        if(params.run_fastp.toString().toBoolean()){
            FASTP_WORKFLOW(sample_ch)
            PHIPPERY(FASTP_WORKFLOW.out)
        
        // Omits Fastp. Assumes fastqs are already trimmed / filtered 
        }else{            
            PHIPPERY(sample_ch)
        }
        PHIPPERYTOAVARDA(PHIPPERY.out)
        virlib = PHIPPERYTOAVARDA.out.virlib
        edgeRhits = PHIPPERYTOAVARDA.out.edgeRhits
        // If params.run_AVARDA was switched on by user, then we run AVARDA 
        // using the outputs from PHIPPERYTOAVARDA. If not, we skip AVARDA 
        if(params.run_AVARDA.toString().toBoolean()){            
            AVARDA(virlib, edgeRhits)
        }
    

    // AVARDA only mode (for VirScan): run AVARDA using user-provided virlib and edgeRhits files, with custom WookScan modifications
    // NO PHIPPERY INVOLVED HERE. USER PROVIDES THEIR OWN VIRLIB AND EDGERHITS, WHICH MAY OR MAY NOT BE DERIVED FROM PHIPPERY OUTPUT
    }else if(params.runtype == 'avarda_only'){
        // We're running AVARDA by itself
        log.info """
        --------------------------------------
        WookScan uses   : AVARDA - AVARDA ONLY
        Runtype         : virscan (but AVARDA only)
        --------------------------------------
        AVARDA is developed by:
        -   Monaco et al.
        Modification performed by:
        -   Preston Leung
        ================================
        virlib          : $params.avarda_names
        edgeRhits       : $params.case_path
        publishDir      : $params.out_path

        """.stripIndent()
        
        virlib = Channel.fromPath(params.avarda_names)        
        edgeRhits =  Channel.fromPath(params.case_path)
        AVARDA(virlib, edgeRhits)
    // HuScan mode: run cutadapt, phippery, and skip AVARDA (incompatible with HuScan), with custom WookScan modifications
    }else if(params.runtype == 'huscan'){
        println "Running HuScan workflow"
        sample_ch = Channel.fromPath(params.sample_table)
        // If we need to process raw paired end files:
        //  -   This will go through CUTADAPT -> Fastp -> PearMerge before going to phippery        
        if(params.run_fastp.toString().toBoolean()){
            println "Running fastp here!!"
            CUTADAPT_WORKFLOW(sample_ch) // Trim the adapters from 5' end in the fastq files            
            FASTP_WORKFLOW(CUTADAPT_WORKFLOW.out) // Filter the fastq files            
            PEARMERGE_WORKFLOW(FASTP_WORKFLOW.out.final_filtered_table_ch) // Merge the paired end reads into single reads
            PHIPPERY(PEARMERGE_WORKFLOW.out.sample_info) // Run Phippery on the merged reads
        }
        // Omits Fastp. Assumes fastqs are already trimmed / filtered (will be added later)
        // }else{            
        //     PHIPPERY(Channel.fromPath(params.sample_table))
        // }
        // // PHIPPERY()        
        // PHIPPERYTOAVARDA(PHIPPERY.out)
        
    }        
    
    // THIS IS WHERE HUSCAN WORKFLOW COULD BRANCH OUT TO DIFFERENT PROCESSES BASED ON FLAGS       
    // SOMEWHERE HERE AT LEAST
    // TODO:
    // We need to figure out if we are performing a VirScan or a HuScan run.
    // Next if it is VirScan, we will proceed as old school. If its HuScan then we need
    // to make it run cutadapt, then fastp, then pear merge and bowtie 2.

  
    

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
