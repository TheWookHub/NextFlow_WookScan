# Nextflow for VirScan Pipeline

**Contributors**
 - Preston
 - Bea
 - Legana
 - Shouyu (Coco) Wei 

### Example command

```
# Running Phippery Only #
nextflow run /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/main.nf \
--run_phippery True \
-profile docker \
--peptide_table InputFiles/peptide_table_VIR3_full_v2.csv \
--sample_table InputFiles/sample_table_UNSW_VirScan.csv \
--run_cpm_enr_workflow true \
--oligo_tile_length 50 \
--read_length 51 \
--results ../TestWookFlow/WookScanNextFlowTest \
--outdir ../TestWookFlow/ \
--dataset_prefix "data"

# Running Phippery & AVARDA #
nextflow run ../nf-core-wookflow/main.nf \
--run_phippery True \
--run_AVARDA True \
-profile docker \
--peptide_table InputFiles/peptide_table_VIR3_full_v2.csv \
--sample_table InputFiles/sample_table_UNSW_VirScan.csv \
--run_cpm_enr_workflow True \
--oligo_tile_length 50 \
--read_length 51 \
--results ../WookScanNextFlowTest2/TestWookFlowPhippery \
--outdir ../WookScanNextFlowTest2/ \
--dict_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/dict_path/blastp_peptide_edges.csv.gz \
--total_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/total_path/sequence_RefSeqOnly_bitscore80plus_total_probability_xr2.csv \
--pairwise_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/pairwise_path/sequence_RefSeqOnly_bitscore80plus_unique_probabilities.csv \
--blast_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/blast_path/sequence_RefSeqOnly_bitscore80plus_bitScore.csv.gz \
--out_path ../WookScanNextFlowTest2/TestWookFlow2_AVARDA \
--out_name PHIPAVARDA

# Running Only AVARDA #
nextflow run ../nf-core-wookflow/main.nf \
--run_AVARDA True \
-profile docker \
--outdir ../WookScanNextFlowTest2/ \
--case ../nf-core-wookflow/subworkflows/local/AVARDA/data/example_input/AVARDA_test_data.tsv.gz \
--dict_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/dict_path/blastp_peptide_edges.csv.gz \
--total_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/total_path/sequence_RefSeqOnly_bitscore80plus_total_probability_xr2.csv \
--pairwise_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/pairwise_path/sequence_RefSeqOnly_bitscore80plus_unique_probabilities.csv \
--blast_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/blast_path/sequence_RefSeqOnly_bitscore80plus_bitScore.csv.gz \
--out_path ../WookScanNextFlowTest2/TestWookFlow2_AVARDA \
--out_name PHIPAVARDA \
--avarda_names ../nf-core-wookflow/subworkflows/local/AVARDA/data/avarda_names/avarda_names.csv.gz

**Note 1**: After gitclone 'https://github.com/TheWookHub/Preston_VirScan_Nextflow.git' in the terminal and entered to **Preston_VirScan_Nextflow** directory. Since the Pre-PhIP-Seq pipeline has not yet been merged into the main branch, the user need to run 
```sh
git branch -a 
```
in the terminal to check which branch you are currently on. Then, run 
```sh
git checkout shouyu_work
```
to move to the branch containing the complete Pre-PhIP-Seq process. After doing this, execute the following command.

**Note 2**: If you want to test the success of the Pre-PhIP-Seq pipeline, you can use the built-in test configuration file. This test configuration file is located in **conf/test.config**. This configuration file uses the virus sequence files located in the **assets/bips_test_data, assets/oligo_csv, assets/oligo_fa** directories of the project. 

Open the **conf/test.config file**, since the test of four modes are commented out, you are required to uncomment them before testing each one. Then, run the pipeline test using the command  
```sh
nextflow run main.nf -profile test
```
This command provides a standardized and fully reproducible method for each workflow branch.

```
# Mode 1: bips_then_dolphyn (Full Workflow) #
# This mode takes a directory of protein sequence files, runs BIPS to generate a barcoded library #
# then uses Dolphyn to generate the predicted epitopes, and pass to BIPS to filter out barcoded library with the predicted epitopes. #
# The input viral sequence file(s) must be .fa or .csv format #
nextflow run /home/shouyu/Preston_VirScan_Nextflow/main.nf \
--mode bips_then_dolphyn \
-profile docker \
--input_viral_seqs_dir /path/to/your/virus_directory/ \
--outdir /path/to/your/result_directory/

# Mode 2: bips_only #
# This mode runs only the BIPS part to generate a complete barcoded oligo library from the input proteomes. #
# The input viral sequence file(s) must be .fa or .csv format #
nextflow run /home/shouyu/Preston_VirScan_Nextflow/main.nf \
--mode bips_only \
-profile docker \
--input_viral_seqs_dir /path/to/your/virus_directory/ \
--outdir /path/to/your/result_directory/

# Mode 3: dolphyn_standalone #
# This mode runs Dolphyn as a standalone tool on protein FASTA file to predict epitopes. #
# The input viral protein sequence file(s) can be .fa single sequences line format or multi sequences lines. dolphyn_standalone can convert multi-line sequences fasta into single line. #
nextflow run /home/shouyu/Preston_VirScan_Nextflow/main.nf \
--mode dolphyn_standalone \
-profile docker  \
--input_protein_fasta_dir /path/to/your/protein_fa_dir \
--outdir /path/to/your/result_directory/


# Mode 4: dolphyn_oligo_only (On Pre-computed Oligos) #
# This mode runs Dolphyn on a set of existing oligonucleotide sequences. #

# Situation A: The input oligos sequence file(s) must be .fa single sequence line format #
nextflow run /home/shouyu/Preston_VirScan_Nextflow/main.nf \
--mode dolphyn_oligo_only \
-profile docker \
--input_oligos_fasta_dir /path/to/your/oligo_fa_dir \
--outdir /path/to/your/result_directory/

# Situation B: The input oligos sequence file(s) must be .csv format #
nextflow run /home/shouyu/Preston_VirScan_Nextflow/main.nf \
--mode dolphyn_oligo_only \
-profile docker \
--input_oligos_csv_dir /path/to/your/oligo_csv_dir \
--outdir /path/to/your/result_directory/

# The Example input oligos sequence files in .fa and .csv are located in "assets/oligo_csv" and "assets/oligo_fa". #

```



## Introduction

**nf-core/wookflow** is a bioinformatics pipeline that integrates phippery and AVARDA together into one nextflow structure. The idea is to allow a smooth flow from fastq files to read counts data that are ready for down stream analyses. While phippery outputs counts data together with edgeR hits to identify peptides that were significantly occurring above background noise, AVARDA will take the hits data to determine which species have been observed based on the peptide hits. AVARDA also accounts for the potential similarity between peptides that come from organisms with high similarity in their genetics. This in turn allows the determination of whether a species can be uniquely identified as due to existence of peptides that were exclusively from that species.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))


> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/wookflow/usage) and the [parameter documentation](https://nf-co.re/wookflow/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/wookflow/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/wookflow/output).

## Credits

nf-core/wookflow was originally written by Preston Leung, but the components that made up the pipeline are written by their respective authors:

1) Phippery - ()
2) AVARDA - ()
3) BuildPhIPSeqLibrary - ()
4) Dolphyn - ()

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#wookflow` channel](https://nfcore.slack.com/channels/wookflow) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/wookflow for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
