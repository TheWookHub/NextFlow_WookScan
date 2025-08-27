# Nextflow for VirScan Pipeline

**Contributors**
 - Preston Leung
 - Bea Delgado-Corrales
 - Legana Fingerhut
 - Shouyu (Coco) Wei 


## Introduction

**nf-core/wookflow** is a bioinformatics pipeline that integrates several tools into one. The pipeline can be broken down into two components. Pre-PhIPSeq oligonucleotide library generation and Post-PhIPSeq sequence analysis. Briefly, PhIPSeq (Phage Immunoprecipitate Sequencing) is a technique that makes use of phages to display peptides on the surface such that antibodies can bind onto. This is followed by separating phages bound by antibodies to those that aren't through the use of magnetic beads. Non-bound beads are then washed away with the remaining bound phages progress to sequencing. For Pre-PhIPSeq BIPS and Dolphyn have been integrated to support custom oligonucleotide library support while phippery supports Post-PhIPSeq data analysis. AVARDA (also Post-PhIPSeq) is a VirScan (defined Human virome oligonucleotide library) Library specific tool that helps to identify individual species of viruses when cross-reactivity exists (Fig. 1).

Pre-PhIPSeq library generation (Fig. 1) aims to provide an integrated approach to go from protein sequences direcctly to oligonucleotide library such that it can be synthesised and be ready for PhIPSeq experiments. Post-PhIPSeq analysis (Fig. 1) aims to allow a smooth flow from fastq files to read counts data that are ready for down stream analyses. While phippery outputs counts data together with edgeR hits to identify peptides that were significantly occurring above background noise, AVARDA will take the hits data to determine which species have been observed based on the peptide hits. AVARDA also accounts for the potential similarity between peptides that come from organisms with high similarity in their genetics. This in turn allows the determination of whether a species can be uniquely identified as due to existence of peptides that were exclusively from that species.


![Test Image](https://github.com/TheWookHub/Preston_VirScan_Nextflow/blob/main_mod/readme_figures/WookFlow_Illustration.png)

***Fig 1.** WookScan overview. Left side illustrates Pre-PhIPSeq library generation and right side shows Post-PhIPSeq analysis. The centre component illustrates the wet laboratory procedure during a PhIPSeq experiment.*

### Example command

#### Running Pre-PhIPSeq pipeline to generate oligo libraries

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

**Note**: If you want to test the success of the Pre-PhIP-Seq pipeline, you can use the built-in test configuration file. This test configuration file is located in **conf/test.config**. This configuration file uses the virus sequence files located in the **assets/bips_test_data, assets/oligo_csv, assets/oligo_fa** directories of the project. 

Open the **conf/test.config file**, since the test of four modes are commented out, you are required to uncomment them before testing each one. Then, run the pipeline test using the command  
```sh
nextflow run main.nf -profile test
```
This command provides a standardized and fully reproducible method for each workflow branch.

#### Running Post-PhIPSeq Analysis on fastq files

General phippery runs when your fastq files are already trimmed and filtered. You can optionally give the raw fastq files to wookflow and then turn `--run_fastp True` in the parameters and Wookflow will run them through fastp. To specify the path of fastq files edit it in the `.csv` file for `--sample_table` input.

```
# Running Phippery Only #

nextflow run /home/preston/PhIPSeq-Pipelines/nf-core-wookflow/main.nf \
--run_phippery True \
-profile docker \
--peptide_table InputFiles/peptide_table_VIR3_full_v7.csv \
--sample_table InputFiles/sample_table_UNSW_VirScan.csv \
--run_cpm_enr_workflow true \
--oligo_tile_length 50 \
--read_length 51 \
--results ../TestWookFlow/WookScanNextFlowTest \
--dataset_prefix "data"

# Running Phippery & AVARDA #

nextflow run ../nf-core-wookflow/main.nf \
--run_phippery True \
--run_AVARDA True \
-profile docker \
--peptide_table InputFiles/peptide_table_VIR3_full_v7.csv \
--sample_table InputFiles/sample_table_UNSW_VirScan.csv \
--run_cpm_enr_workflow True \
--oligo_tile_length 50 \
--read_length 51 \
--results ../WookScanNextFlowTest2/phippery_results \
--out_path ../WookScanNextFlowTest2/avarda_reults \
--out_name PHIPAVARDA \
--max_cpus 24 \
--max_memory '36.GB' \
--run_fastp True


# Running Only AVARDA (assumes you have some hits data) #
# Note: AVARDA only works for VirScan Library for now.

nextflow run ../nf-core-wookflow/main.nf \
--run_AVARDA True \
-profile docker \
--case_path ../nf-core-wookflow/subworkflows/local/AVARDA/data/example_input/AVARDA_test_data.tsv.gz \
--avarda_names ../nf-core-wookflow/subworkflows/local/AVARDA/data/avarda_names/avarda_names.csv.gz
--out_path ../WookScanNextFlowTest2/TestWookFlow2_AVARDA \
--out_name PHIPAVARDA \
--max_cpus 24 \
--max_memory '36.GB'

``` 



## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/wookflow/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/wookflow/output).

## Credits

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#wookflow` channel](https://nfcore.slack.com/channels/wookflow) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/wookflow for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

nf-core/wookflow was originally written by Preston Leung, but the components that made up the pipeline are written by their respective authors. An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.




## Other stuff todo later

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/wookflow/usage) and the [parameter documentation](https://nf-co.re/wookflow/parameters).



