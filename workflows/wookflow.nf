/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { VIRALDB                } from '../subworkflows/local/viraldb/main.nf'
include { FASTP_WORKFLOW          } from '../subworkflows/local/fastp/main.nf'
include { PHIPPERY                } from '../subworkflows/local/phippery/main.nf'
include { PHIPPERYTOAVARDA        } from '../subworkflows/local/phipperyToAvarda/main.nf'
include { AVARDA                  } from '../subworkflows/local/AVARDA/main.nf'

// may remove this later
// include { POSTAVARDA_WORKFLOW     } from '../subworkflows/local/PostAVARDA/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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


workflow WOOKFLOW {    
    main:   

    //
    // RUN: ViralDB - Build Viral Library support files for AVARDA        
    // TODO   

    
    // Print out the message for running phippert    
    if(params.run_phippery){        
        log.info """
            --------------------------------------
            WookScan uses: P H I P - F L O W!
            --------------------------------------            
            Phippery & phip-flow is developed by:
            -   Matsen, Overbaugh, and Minot Labs
            -   Fred Hutchinson CRC, Seattle WA
            Repository: 
            -   https://github.com/matsengrp/phip-flow            
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
            -   Preston Leung
            ================================
            user_pep_id     : $params.user_pep_id
            publishDir      : $params.results

        """.stripIndent()        
    }
    
    // Print out message for AVARDA with/without PHIPPERY messages    
    if(params.run_AVARDA){        
        if(params.run_phippery){
            // If we take stuff directly from phippery output
            log.info """
            --------------------------------------
            WookScan uses: AVARDA: PHIPPERY-AVARDA
            --------------------------------------
            AVARDA is developed by:
            -   Monaco et al.
            Repository: 
            -   https://github.com/drmonaco/AVARDA
            Modification performed by:
            -   Preston Leung
            ================================
            virlib      : From PHIPPERYTOAVARDA.out.virlib
            edgeRhits   : From PHIPPERYTOAVARDA.out.edgeRhits
            publishDir  : $params.out_path
            
            """.stripIndent()            
        }else{
            // We're running AVARDA by itself
            log.info """
            --------------------------------------
            WookScan uses: AVARDA: AVARDA ONLY
            --------------------------------------
            AVARDA is developed by:
            -   Monaco et al.
            Modification performed by:
            -   Preston Leung
            ================================
            virlib      : $params.avarda_names
            edgeRhits   : $params.case_path
            publishDir  : $params.out_path

            """.stripIndent()
            
        }
    }
    // Running Fastp -> Phippery
    if(params.run_phippery){
        FASTP_WORKFLOW()
        PHIPPERY(FASTP_WORKFLOW.out)
        // PHIPPERY()        
        PHIPPERYTOAVARDA(PHIPPERY.out)
    }
    // Running AVARDA
    if(params.run_AVARDA){ 
        if(params.run_phippery){
            virlib = PHIPPERYTOAVARDA.out.virlib
            edgeRhits = PHIPPERYTOAVARDA.out.edgeRhits
        }else{
            virlib = Channel.fromPath(params.avarda_names)
            edgeRhits =  Channel.fromPath(params.case_path)        
        }        
        AVARDA(virlib, edgeRhits)
        // probably don't need this anymore since the sample names are correclty labeled in the AVARDA output
        // POSTAVARDA_WORKFLOW() 
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
