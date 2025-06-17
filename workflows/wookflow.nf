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
include { FASTP_WORKFLOW          } from '../subworkflows/local/fastp/main.nf'
include { PHIPPERY                } from '../subworkflows/local/phippery/main.nf'
include { PHIPPERYTOAVARDA        } from '../subworkflows/local/phipperyToAvarda/main.nf'
include { AVARDA                  } from '../subworkflows/local/AVARDA/main.nf'
include { POSTAVARDA_WORKFLOW     } from '../subworkflows/local/PostAVARDA/main.nf'
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
    // phippery can read the trimmed fastq output from fastp. Also figure out where
    // to put this part in the pipeline.


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
        FASTP_WORKFLOW()
        // PHIPPERY()
        // PHIPPERYTOAVARDA(PHIPPERY.out)
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
        // AVARDA(virlib, edgeRhits)
        // POSTAVARDA_WORKFLOW()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
