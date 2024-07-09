#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

formatted_sumstats <- as.character(args[1])
gwas_name <- as.character(args[2])
tissue <- as.character(args[3])
total_chunks <- as.numeric(args[4])
chunk <- as.numeric(args[5])
fusion_weight_dir  <- as.character(args[6])

library(RSQLite)
library(ctwas)
library(data.table)

z_snp <- readRDS(formatted_sumstats)

#ld_pgenfs <- paste0("~/cTWAS/data/GTEx_LDref/", paste0("GTEx_EUR_chr", 1:22, ".bed"))
##ld_R_dir <- "/home/jupyter/cTWAS/data/GTEx_LDref/matrices_0.1/LDR_b38_cova/"
ld_R_dir <- "/home/jupyter/cTWAS/data/UKB_LDref/matrices_0.1/LDR_b38_cova/"

outputdir <- paste0("~/cTWAS/results/imputed_gene_expression/", gwas_name, "/", tissue, "/chunk", chunk)
outname <- paste0(gwas_name, "_chunk", chunk, "_UKBld")
dir.create(outputdir, showWarnings=F, recursive=T)

fusion_weight_dir_chunk <- paste0(fusion_weight_dir, '/chunk', chunk, '/')

# get gene z score
res <- impute_expr_z(z_snp = z_snp,
                     weight = fusion_weight_dir_chunk,
                     method = "blup", 
                     ld_R_dir = ld_R_dir,
                     outputdir = outputdir,
                     outname = outname,
                     harmonize_z = T,
                     harmonize_wgt = T,
                     strand_ambig_action_z = "drop",
                     recover_strand_ambig_wgt = T)

z_snp <- res$z_snp
z_gene <- res$z_gene

save(z_snp, file=paste0(outputdir, '/', gwas_name, '_chunk', chunk, '_UKBld_Zsnp_results.Rd'))
save(z_gene, file=paste0(outputdir, '/', gwas_name, '_chunk', chunk, '_UKBld_Zgene_results.Rd'))
