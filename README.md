# Nextflow for VirScan Pipeline

TO BE FILLED!

**Contributors**
 - Preston
 - Bea
 - Legana


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
```



## Introduction

**nf-core/wookflow** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/wookflow \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/wookflow/usage) and the [parameter documentation](https://nf-co.re/wookflow/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/wookflow/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/wookflow/output).

## Credits

nf-core/wookflow was originally written by Preston Leung.

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
