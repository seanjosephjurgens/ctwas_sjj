#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

base_wd <- "/home/jupyter/cTWAS/"

formatted_sumstats <- as.character(args[1])
gwas_name <- as.character(args[2])
tissue <- as.character(args[3])
total_chunks <- as.numeric(args[4])
chunk <- as.numeric(args[5])
try(base_wd <- as.character(args[6]))

library(RSQLite)
library(ctwas)
library(data.table)

z_snp <- readRDS(formatted_sumstats)

#### formatting prediction models #### 
#specify the weight to remove lncRNA from
weight <- paste0(base_wd, "data/GTEx_models/weights/eqtl/mashr/mashr_", tissue, ".db")
#read the PredictDB weights
sqlite <- dbDriver("SQLite")
db = dbConnect(sqlite, weight)
query <- function(...) dbGetQuery(db, ...)
weights_table <- query("select * from weights")
extra_table <- query("select * from extra")
dbDisconnect(db)
#subset to protein coding genes only
extra_table <-  extra_table[extra_table$gene_type=="protein_coding",,drop=F]
weights_table <- weights_table[weights_table$gene %in% extra_table$gene,]
#subset to chunk only
genes <- extra_table$gene
genes_split <- split(genes,  cut(seq_along(genes), total_chunks, labels = FALSE))
gene_subset <- genes_split[[chunk]]
extra_table <-  extra_table[extra_table$gene %in% gene_subset,,drop=F]
weights_table <- weights_table[weights_table$gene %in% extra_table$gene,]
#subset the covariances
weight_info = read.table(gzfile(paste0(tools::file_path_sans_ext(weight), ".txt.gz")), header = T)
weight_info <- weight_info[weight_info$GENE %in% extra_table$gene,]
#write the .db file and the covariances
if (!file.exists(paste0(base_wd, "data/GTEx_models/weights_nolnc/mashr_", tissue, "_nolnc_subset_chunk", chunk, ".db"))){
  db <- dbConnect(sqlite, paste0(base_wd, "data/GTEx_models/weights_nolnc/mashr_", tissue, "_nolnc_subset_chunk", chunk, ".db"))
  dbWriteTable(db, "extra", extra_table)
  dbWriteTable(db, "weights", weights_table)
  dbDisconnect(db)
  
  weight_info_gz <- gzfile(paste0(base_wd, "data/GTEx_models/weights_nolnc/mashr_", tissue, "_nolnc_subset_chunk", chunk, ".txt.gz"), "w")
  write.table(weight_info, weight_info_gz, sep=" ", quote=F, row.names=F, col.names=T)
  close(weight_info_gz)
}
#specify the weight for the analysis
weight_subset <- paste0(base_wd, "data/GTEx_models/weights_nolnc/mashr_", tissue, "_nolnc_subset_chunk", chunk, ".db")

#ld_pgenfs <- paste0("~/cTWAS/data/GTEx_LDref/", paste0("GTEx_EUR_chr", 1:22, ".bed"))
ld_R_dir <- paste0(base_wd, "data/GTEx_LDref/matrices_0.1/LDR_b38_cova/")

outputdir <- paste0(base_wd, "results/imputed_gene_expression/", gwas_name, "/", tissue, "/chunk", chunk)
outname <- paste0(gwas_name, "_chunk", chunk)
dir.create(outputdir, showWarnings=F, recursive=T)

# get gene z score
res <- impute_expr_z(z_snp = z_snp,
                     weight = weight_subset,
                     ld_R_dir = ld_R_dir,
                     outputdir = outputdir,
                     outname = outname,
                     harmonize_z = T,
                     harmonize_wgt = T,
                     strand_ambig_action_z = "none",
                     recover_strand_ambig_wgt = T)

z_snp <- res$z_snp,
save(z_snp, file=paste0(outputdir, '/', gwas_name, '_chunk', chunk, '_Zsnp_results.Rd'))
z_gene <- res$z_gene,
save(z_gene, file=paste0(outputdir, '/', gwas_name, '_chunk', chunk, '_Zgene_results.Rd'))
###save(res, file=paste0(outputdir, '/', gwas_name, '_chunk', chunk, '_Zresults.Rd'))

system(paste0("rm -rf ", base_wd, "data/GTEx_models/weights_nolnc/mashr_", tissue, "_nolnc_subset_chunk", chunk, ".db"))
system(paste0("rm -rf ", base_wd, "data/GTEx_models/weights_nolnc/mashr_", tissue, "_nolnc_subset_chunk", chunk, ".txt.gz"))
