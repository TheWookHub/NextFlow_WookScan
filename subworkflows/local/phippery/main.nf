/*
 * This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENCE
 * License, v. 3.0. 
 */


/* 
 * 'PhIP-Flow' - A Nextflow pipeline for running common phip-seq analysis workflows
 * 
 * Fred Hutchinson Cancer Research Center, Seattle WA.
 * 
 * Jared Galloway
 * Kevin Sung
 * Sam Minot
 * Erick Matsen 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters - example data get's run by default
 # $baseDir 
 */ 
// params.sample_table     = "$baseDir/data/pan-cov-example/sample_table_with_beads_and_lib.csv"
// if (params.sample_table != "$baseDir/data/pan-cov-example/sample_table_with_beads_and_lib.csv")
//     params.reads_prefix = "$launchDir"
// else
//     params.reads_prefix = "$baseDir"
// params.peptide_table    = "$baseDir/data/pan-cov-example/peptide_table.csv"
// params.results          = "$PWD/results/"


// params.sample_table     = "./data/pan-cov-example/sample_table_with_beads_and_lib.csv"
// if (params.sample_table != "./data/pan-cov-example/sample_table_with_beads_and_lib.csv")
//     params.reads_prefix = "$launchDir"
// else
//     params.reads_prefix = "$baseDir"
// params.peptide_table    = "./data/pan-cov-example/peptide_table.csv"
// params.results          = "$PWD/results/"


// log.info """\
// P H I P - F L O W!
// Matsen, Overbaugh, and Minot Labs
// Fred Hutchinson CRC, Seattle WA
// ================================
// sample_table    : $params.sample_table
// peptide_table   : $params.peptide_table
// results         : $params.results
// reads_prefix    : $params.reads_prefix

// """

/* 
 * Import modules 
 */
nextflow.enable.dsl=2

include { ALIGN } from './workflows/alignment.nf'
include { STATS } from './workflows/statistics.nf'
include { DSOUT } from './workflows/output.nf'
include { AGG } from './workflows/aggregate.nf'

workflow PHIPPERY {    
    take:
        filtered_sample_table_ch    
    
    main:
    ALIGN(filtered_sample_table_ch) | STATS | DSOUT | AGG     
        // original flow
    // ALIGN | STATS | DSOUT | AGG

    emit:
    AGG.out
}
