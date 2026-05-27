#!/bin/bash

set -euo pipefail

FASTA=!{oligo_fasta}
CPUS=!{task.cpus}

mkdir peptide_index
bowtie2-build \
    --threads $CPUS \
    $FASTA \
    peptide_index/peptide
