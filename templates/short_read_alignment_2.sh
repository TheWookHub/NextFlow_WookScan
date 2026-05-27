#!/bin/bash

: '
Alignment of oligonucleotides to reference library using bowtie2.

Parameters set for bowtie2 are as close to bowtie as possible.

We report only the best alignment found, and no more.
'

set -euo pipefail

STREAM_FILE_CMD=!{params.fastq_stream_func}
FASTQ=!{respective_replicate_path}
INDEX=!{index}/peptide
ALIGN_OUT_FN=!{sample_id}.sam
# READ_LENGTH=!{params.read_length}
# PEPTIDE_LENGTH=!{params.oligo_tile_length}
CPUS=!{task.cpus}
MM=!{params.n_mismatches}
OP_ARGS="!{params.bowtie_2_optional_args}" # because this is a whole bunch of optional args so need "" to make it a single string.
TRIM3=!{params.bowtie_2_trim3}

# if [ ${PEPTIDE_LENGTH} -lt ${READ_LENGTH} ]; then
#     let TRIM3=${READ_LENGTH}-${PEPTIDE_LENGTH}
# else
#     TRIM3=0
# fi

echo $OP_ARGS

$STREAM_FILE_CMD $FASTQ | bowtie2 \
  --trim3 $TRIM3 \
  --threads $CPUS \
  -x $INDEX \
  $OP_ARGS \
  -S $ALIGN_OUT_FN \
  -U -
  
