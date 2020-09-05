#!/usr/bin/env Rscript
library(sleuth)
library(tidyverse)
library(biomaRt)

args <- commandArgs(TRUE)

sample_id <- dir(args[1])
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(args[1], id))
kal_dirs

s2c <- readr::read_tsv(args[2], col_names = TRUE) 
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)


# Extract gene or transcripts infromation from Ensemble data source
ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Load the kallisto processed data into the object
so_geneMode <- sleuth_prep(s2c, 
	target_mapping = t2g, 
	read_bootstrap_tpm = TRUE, 
	extra_bootstrap_summary = TRUE, 
	aggregation_column = 'ens_gene',
	gene_mode = TRUE
)

# fit models for significant differentially expressed genes
so_geneMode <- sleuth_fit(so_geneMode, ~condition, 'full')
so_geneMode <- sleuth_fit(so_geneMode, ~1, 'reduced')
so_geneMode <- sleuth_wt(so_geneMode, 'conditioncontrol')
so_geneMode <- sleuth_lrt(so_geneMode, 'reduced', 'full')

# export result table
full_results <- sleuth_results(so_geneMode, 'reduced:full', 'lrt', show_all = FALSE)
readr::write_tsv(full_results, "full_sleuth_results.tsv", col_names = TRUE)

sleuth_significant <- dplyr::filter(full_results, qval <= 0.05)
readr::write_tsv(sleuth_significant, "sleuth_significant.tsv", col_names = TRUE)

# Generate the whole gene table

gene_info <- full_results %>%
	dplyr::select(target_id, ext_gene)

gene_table <- kallisto_table(so_geneMode) %>%
	dplyr::left_join(., gene_info, by = "target_id") %>%
	dplyr::select(sample, target_id, ext_gene, everything())

readr::write_tsv(gene_table, "gene_table.tsv", col_names = TRUE)
# Saving for shinny
save(so_geneMode, file=paste("geneMode.so"))
