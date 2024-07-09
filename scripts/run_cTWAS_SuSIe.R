#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

base_wd <- "/home/jupyter/cTWAS/"

gwas_name <- as.character(args[1])
pheno_name <- as.character(args[2])
p_thres <- as.numeric(args[3])
r2_thres <- as.numeric(args[4])
rerun <- as.logical(args[5])
try(base_wd <- as.character(args[6]))

####################
## Collect results
####################

library(data.table)

#gwas_name='DCM_GWAS'
#pheno_name='GTEx_Heart_Left_Ventricle'
#p_thres=p1e-07
#r2_thres=r20.001

finemap_dir <- paste0(base_wd, 'results/fine_mapping/')
tissue <- paste0(pheno_name, "_p", p_thres, "_r2", r2_thres)
outputdir <- paste0(finemap_dir, '/', gwas_name, '/', tissue, '/')
outname <- paste0(gwas_name, '_', tissue, '_finemap_results_UKBld')
res_file <- paste0(outputdir, "/", outname)

#message(res_file)
if( !file.exists(paste0(res_file, ".susieIrss.txt")) ){
if( !rerun ){

new_dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/combined/')
system(paste0("mkdir ", new_dir))

# Collect per-chromosome results for each chunk and munge
for(chr in c(1:22)){
  
  message("Busy with chr", chr)
  
  qclist_total<- list()
  wgtlist_total <- list()
  z_gene_chr_total <- NULL
  exprvar_total <- NULL
  
  for(chunk in c(1:65)){
    dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/chunk', chunk)
    setwd(dir)
    expr_file <- paste0(gwas_name, '_chunk', chunk, '_UKBld_chr', chr, '.exprqc.Rd')
    load(expr_file)
    qclist_total <- c(qclist_total, qclist)
    wgtlist_total <- c(wgtlist_total, wgtlist)
    z_gene_chr_total <- rbind(z_gene_chr_total, z_gene_chr)
    
    exprvar_file <- paste0(gwas_name, '_chunk', chunk, '_UKBld_chr', chr, '.exprvar')
    exprvar <- fread(exprvar_file, stringsAsFactors = F, data.table = F)
    exprvar_total <- rbind(exprvar_total, exprvar)
  }
  
  qclist <- qclist_total 
  wgtlist <- wgtlist_total
  z_gene_chr <- z_gene_chr_total
  
  fail1 <- length(which(is.na(z_gene_chr$id)))
  message("\t", fail1, " number of NA gene ids")
  fail2 <- length(which(is.na(z_gene_chr$z)))
  message("\t", fail2, " number of NA gene imputations")
  
  if(fail1>0 | fail2>0){
    message("removing failed genes")
    if(fail1>0){
      z_gene_chr <- z_gene_chr[-which(is.na(z_gene_chr$id)), ]
      fail2 <- length(which(is.na(z_gene_chr$z)))
    }
    if(fail2>0){
      failed_genes <- z_gene_chr[which(is.na(z_gene_chr$z)), 'id']
      z_gene_chr <- z_gene_chr[-which(z_gene_chr$id %in% failed_genes), ]
      qclist <- qclist[-which(names(qclist) %in% failed_genes)]
      wgtlist <- wgtlist[-which(names(wgtlist) %in% failed_genes)]
      exprvar_total <- exprvar_total[-which(exprvar_total$id %in% failed_genes), ]
    }
    
  }
  
  fail <- FALSE
  for(i in c(1:length(qclist))){
    fail <- c(fail, any(is.na(qclist[[i]])))
  }
  if(any(fail)){
    message("\tWARNING some NA values in qclist file")
  }
  
  fail <- FALSE
  for(i in c(1:length(wgtlist))){
    fail <- c(fail, any(is.na(wgtlist[[i]])))
  }
  if(any(fail)){
    message("\tWARNING some NA values in wgtlist file")
  }
  
  if(any(is.na(exprvar_total))){
    message("\tWARNING some NA values in exprvar_total file")
  }
  
  outfile1 <- paste0(gwas_name, '_chunkall_UKBld_chr', chr, '.exprqc.Rd')
  save(qclist, wgtlist, z_gene_chr, 
       file=paste0(new_dir, '/', outfile1)
  )
  outfile2 <- paste0(gwas_name, '_chunkall_UKBld_chr', chr, '.exprvar')
  write.table(exprvar_total, file=paste0(new_dir, '/', outfile2), col.names=T, row.names=F, quote=F, sep='\t')
  
}

# Collect LD regions (same for every chunk, so just use chrunk1 results)
for(chr in c(1:22)){
  dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/chunk1/')
  new_dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/combined/')
  
  old_file <- paste0(gwas_name, '_chunk1_UKBld_ld_R_chr', chr, '.txt')
  new_file <- paste0(gwas_name, '_chunkall_UKBld_ld_R_chr', chr, '.txt')
  
  system(paste0('cp ', dir, '/', old_file, ' ', new_dir, '/', new_file))
}

# Collect main Zscore results (one result per junk)
chunk <- 1
dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/chunk', chunk)
setwd(dir)
zfile <- paste0(gwas_name, '_chunk', chunk, '_UKBld_Zsnp_results.Rd')
load(zfile)
z_snp_total <- z_snp

z_gene_total <- NULL
for(chunk in c(1:65)){
  message(paste0("Busy with ", chunk, "..."))
  dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/chunk', chunk)
  setwd(dir)
  zfile <- paste0(gwas_name, '_chunk', chunk, '_UKBld_Zgene_results.Rd')
  load(zfile)
  z_gene_total <- rbind(z_gene_total, z_gene)
}
z_snp <- z_snp_total
z_gene <- z_gene_total

rm <- which(is.na(z_snp$z))
message(length(rm), " missing SNP Z scores.")
if(length(rm)>0){
  message('removing those...')
  z_snp <- z_snp[-rm, ]
}
rm <- which(is.na(z_gene$z))
message(length(rm), " missing GENE Z scores.")
if(length(rm)>0){
  message('removing those...')
  z_gene <- z_gene[-rm, ]
}

new_dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/combined/')
setwd(new_dir)
res <- list(z_gene, z_snp)
save(res, file=paste0(gwas_name, '_chunkall_UKBld_Zresults.Rd'))


################################
#### Remove some unneeded files 
################################

z_gene_total <- NULL
for(chunk in c(1:65)){
  message(paste0("Busy with ", chunk, "..."))
  dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/chunk', chunk)
  setwd(dir)
  zfile <- paste0(gwas_name, '_chunk', chunk, '_UKBld_Zgene_results.Rd')
  system(paste0("rm  ", zfile))
  zfile <- paste0(gwas_name, '_chunk', chunk, '_UKBld_Zsnp_results.Rd')
  system(paste0("rm  ", zfile))
}

}

################################
# Run cTWAS: FineMapping
################################

library(ctwas)
library(data.table)
new_dir <- paste0(base_wd, 'results/imputed_gene_expression/', gwas_name, '/', pheno_name, '_p', p_thres, '_r2', r2_thres, '/combined/')
for(chr in c(1:22)){
  outfile1 <- paste0(gwas_name, '_chunkall_UKBld_chr', chr, '.exprqc.Rd')
  load(paste0(new_dir, "/", outfile1))
}

res <- get(load(paste0(new_dir, '/', gwas_name, '_chunkall_UKBld_Zresults.Rd')))
z_gene <- res[[1]]
z_snp <- res[[2]]

ld_exprvarfs <- paste0(new_dir, '/', gwas_name, '_chunkall_UKBld_chr', c(1:22), '.exprvar')
ld_R_dir <- paste0(base_wd, "data/UKB_LDref/matrices_0.1/LDR_b38_cova/")
#ld_pgenfs <- paste0("/home/jupyter/cTWAS/data/GTEx_LDref/GTEx_EUR_chr", 1:22, ".bed")

regions <- system.file("extdata/ldetect", "EUR.b38.bed", package = "ctwas")
regions_df <- read.table(regions, header = T)

tissue <- paste0(pheno_name, '_p', p_thres, '_r2', r2_thres)
finemap_dir <- paste0(base_wd, 'results/fine_mapping/')
system(paste0("mkdir ", finemap_dir))
system(paste0("mkdir ", finemap_dir, '/', gwas_name, '/'))
system(paste0("mkdir ", finemap_dir, '/', gwas_name, '/', tissue, '/'))
outputdir <- paste0(finemap_dir, '/', gwas_name, '/', tissue, '/')
outname <- paste0(gwas_name, '_', tissue, '_finemap_results_UKBld')

# run ctwas_rss
detach("package:ctwas", unload=TRUE)
remotes::install_github("seanjosephjurgens/ctwas_sjj",ref = "main")

library(ctwas)
ncore <- 16
ncore.rerun <- 4
ncore_index_regions <- 11
niter1 <- 3
niter2 <- 30
thin <- 0.1
max_snp_region <- 20000
ctwas_rss(z_gene = z_gene, 
          z_snp = z_snp, 
          ld_exprvarfs = ld_exprvarfs, 
          ld_R_dir = ld_R_dir, 
          #ld_pgenfs = ld_pgenfs,
          ld_regions = "EUR",
          ld_regions_version = "b38",
          outputdir = outputdir, 
          outname = outname,
          thin = thin,
          max_snp_region = max_snp_region,
          ncore = ncore,
          ncore.rerun = ncore.rerun,
          ncore_index_regions = ncore_index_regions,
          
          prob_single = 0.5,
          rerun_gene_PIP = 0.5,
          niter1 = niter1,
          niter2 = niter2,
          L = 5,
          group_prior = NULL,
          group_prior_var = NULL,
          estimate_group_prior = T,
          estimate_group_prior_var = T,
)

system(paste0("rm -r ", outputdir, outname, "_LDR"))

}else{
	
	message("output file:")
	message(paste0("'", res_file, "'"))
	message("output file exists. no rerun needed. stopping.")

}
