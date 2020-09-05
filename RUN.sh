#!/usr/bin/env bash

nextflow run main.nf -profile standard -resume --fqs="/home/group2/LPL/pipeline/RNAseq-nf/sample_sheet.tsv" --experiment="/home/group2/LPL/pipeline/RNAseq-nf/info.txt" --transcriptome="/home/group2/LPL/reference_database/hg38/transcriptome/transcriptome.fa"
