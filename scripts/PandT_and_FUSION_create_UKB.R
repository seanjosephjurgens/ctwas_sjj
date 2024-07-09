#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)

########################
#Prune-and-thresholding
########################

####### P+T in UKBB using PLINK

#p_thresholds <- c(1e-7)
#r2_thresholds <- c(0.0005)
#pheno_name <- 'eQTLGen'
#base_wd <- '~/cTWAS/data/eQTLGen'
#pheno_file <- paste0('~/cTWAS/data/eQTLGen/eQTLGen_eQTLs_withbetasandses_atleast40percsamplesize_liftedB38_QCed.tsv.gz')
#cis_trans <- "cis"
#gwas_file <- paste0('~/cTWAS/data/sum_stats/cTWAS_input/DCM_GWAS_nonimputed_SS.RDS')
#gwas_name <- "DCM_GWAS"
#run_MR <- TRUE

p_thresholds <- c(as.numeric(args[1]))
r2_thresholds <- c(as.numeric(args[2]))
pheno_name <- as.character(args[3])
base_wd <- as.character(args[4])
pheno_file <- as.character(args[5])
cis_trans <- as.character(args[6])
gwas_file <- as.character(args[7])
gwas_name <- as.character(args[8])
run_MR <- as.logical(args[9])

gwas_file <- gsub("SS", "SS_BETA_SE", gwas_file)

pheno <- paste0(gwas_name, "__", pheno_name)
wd <- paste0(base_wd, '/', pheno)
try(system(paste0("mkdir ", wd)))

message('\npheno is ',pheno,'.\n')
message(paste0('reading GWAS data from: ', gwas_file, '\n'))

tot <- fread(paste0(pheno_file), stringsAsFactors=F)
if(cis_trans=="cis" & "cis_trans" %in% colnames(tot)){
  tot <- tot[tot$cis_trans=="cis", ]
}
tot <- as.data.frame(tot)
tot$CHROM <- as.numeric(tot$CHROM)

wd <- paste0(wd,'/P_T/')
try(system(paste0("mkdir ", wd)))
setwd(wd)

### Filter to variants found in the respective GWAS dataset
gwas <- readRDS(paste0(gwas_file))
tot <- tot[tot$SNP %in% gwas$id, ]

library(doParallel)
library(parallel)
cl <- parallel::makeForkCluster(22, outfile = "")
doParallel::registerDoParallel(cl)

foreach (CHR = c(1:22), .combine='c',
         .packages = "ctwas") %dopar% {
           
  cat("Busy with chromosome ", CHR, "....\n")
  cat("\n\n\n\n\n")
  Sys.sleep(4)
      
  LDreference_prefix <- paste0(base_wd, '/../UKB_LDref/UKBB_ldref_chr',CHR,'_sorted')
  library(data.table)
  
  
  dat <- tot[tot$CHROM==CHR,]
  if(nrow(dat)>0){
    
    # Making Clump readable file (we will use LDpred for P+T method)
    clump_results <- list(NULL)
    cat('Starting Clumping process..\n')
    ntotal <- length(unique(dat$ensembl_id))
    start <- 1
    for(i in unique(dat$ensembl_id)){
      cat('Busy with gene', i, 'which is gene number', start, 'of',ntotal,'genes...\n')
      start <- start + 1
      Clumpfile <- dat[dat$ensembl_id==i,]
      write.table(Clumpfile[,c("SNP", "ALLELE1", "slope", "pval_nominal")], file=paste0(wd,pheno,'_Clumpfile_', i,'.txt'), col.names=T, row.names=F, quote=F)
      
      for(r2 in r2_thresholds){
        for(p in p_thresholds){
          cat('   ... p==',p,'...\n')
          system(paste0('~/plink2 --pfile ',LDreference_prefix,
                        ' --clump ',wd,pheno,'_Clumpfile_', i,'.txt --clump-unphased --clump-field pval_nominal --clump-p1 ', p, ' --clump-p2 ', p,
                        ' --clump-r2  ', r2, '  --clump-kb 10000 --out ',wd,pheno,'_',i,'_Clumped_p',p,'_r2', r2, '.txt'), intern=T)
          if(file.exists(paste0(wd,pheno,'_',i,'_Clumped_p',p,'_r2', r2, '.txt.clumps'))){
            interres <- fread(paste0(wd,pheno,'_',i,'_Clumped_p',p,'_r2', r2, '.txt.clumps'), stringsAsFactors=F)
            interres <- as.data.frame(interres$ID)
            colnames(interres)<-"SNP"
            Clumpedfile <- merge(interres, Clumpfile, by="SNP", all=F)
            system(paste0('rm ', wd,pheno,'_',i,'_Clumped_p',p,'_r2', r2, '.txt.clump*'))
            clump_results[[paste0('p_',p,'_r2_', r2)]] <- rbind(clump_results[[paste0('p_',p,'_r2_',r2)]], Clumpedfile)
          }
          system(paste0('rm ', wd,pheno,'_',i,'_Clumped_p',p,'_r2', r2, '.txt.log'))
          message(paste0(round(start/ntotal*100,2)), " completed...")
        }}
      system(paste0('rm ',wd,pheno,'_Clumpfile_', i,'.txt'))
    }
    
    for(r2 in r2_thresholds){
      for(p in p_thresholds){
        if(nrow(clump_results[[paste0('p_',p,'_r2_', r2)]])>0){
          write.table(clump_results[[paste0('p_',p,'_r2_', r2)]], file=paste0(wd,pheno,'_allgenes_Clumped10Mb_chr',CHR,'_p',p,'_r2', r2, '.txt'), col.names=T, row.names=F, quote=F, sep='\t')
        }
      }}
  }else{
    cat('No variant-gene associations found for chr', CHR, '!\n')
  }

}
parallel::stopCluster(cl)

########################
# Make FUSION pred mods
########################

### Example:
#load('/home/jupyter/packages/ctwas/extdata/example_fusion_weights/Tissue/Tissue.gene1.wgt.RDat')
#pos <- fread('/home/jupyter/packages/ctwas/extdata/example_fusion_weights/Tissue.pos')

wd2 <- sub("P_T", "FUSION_weights", wd, )
system(paste0("mkdir ", wd2))

for(r2 in r2_thresholds){
  for(p in p_thresholds){

    wd3 <- paste0(wd2, "/p",p,'_r2', r2, '/')
    system(paste0("mkdir ", wd3))
        
    clumped_all <- NULL
    for(CHR in c(1:22)){
      clumped <- fread(paste0(wd, pheno,'_allgenes_Clumped10Mb_chr',CHR,'_p',p,'_r2', r2, '.txt'), stringsAsFactors = F, data.table=F)
      clumped_all <- rbind(clumped_all, clumped)                
    }
    ## Make 65 chunks
    genes <- unique(clumped_all$ensembl_id)
    total_chunks <- 65
    genes_split <- split(genes,  cut(seq_along(genes), total_chunks, labels = FALSE))
    
    cl <- parallel::makeForkCluster(65, outfile = "")
    doParallel::registerDoParallel(cl)
    
    foreach (j_num = c(1:65), .combine='c', .packages = "ctwas") %dopar% {
               
      j <- genes_split[[j_num]]

      wd4 <- paste0(wd3, "/chunk", j_num, "/")
      system(paste0("mkdir ", wd4), ignore.stderr=T)
      
      j <- as.vector(unlist(j))
      pos_file <- NULL
	
      MR_exposures <- NULL
      for(i in j){
        #cat("Busy with ", i, "...\n")
        clumped <- clumped_all[clumped_all$ensembl_id==i, ]
        clumped <- clumped[order(clumped$pval_nominal), ]
        if(length(which(duplicated(clumped$SNP)))>0){
          clumped <- clumped[-(which(duplicated(clumped$SNP))), ]
        }
        clumped <- clumped[order(clumped$GENPOS), ]
        rownames(clumped) <- clumped$SNP
        if(!"Z" %in% colnames(clumped)){
          clumped$Z <- clumped$slope / clumped$slope_se
        }
        wgt.matrix <- clumped[,c("Z", "Z", "Z", "Z")]
        colnames(wgt.matrix) <- c("blup", "lasso", "top1", "enet")
        wgt.matrix <- as.matrix(wgt.matrix)
        if(nrow(clumped)>1){
          clumped$order <- c(1:nrow(clumped))
        }else{
          clumped$order <- 1
        }  
        rownames(clumped) <- clumped$order
        clumped$empty <- 0
        snps <- clumped[, c("CHROM", "SNP","empty",  "GENPOS", "ALLELE1", "ALLELE0")]
        colnames(snps) <- c("V1", "V2", "V3", "V4", "V5", "V6")
        outfile <- paste0(i, '.wgt.RDat')
        
        if(any(is.na(wgt.matrix)) | any(is.na(snps))){
          message("potential issue for gene", i, "!")
          print(wgt.matrix)
          print(snps)
        }
        
        if(run_MR){
          MR_clumped <- clumped[,c('ensembl_id', 'SNP', 'CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'slope', 'slope_se', 'pval_nominal', 'N')]
          MR_exposures <- rbind(MR_exposures, MR_clumped)
        }
        save(wgt.matrix, snps, file=paste0(wd4, outfile))
        
        pos_line <- c(pheno, paste0("chunk", j_num, '/', outfile), i, clumped[1,'CHROM'], min(clumped$GENPOS)-1000, max(clumped$GENPOS)+1000, max(clumped$N))
        pos_file <- rbind(pos_file, pos_line)
      }
      colnames(pos_file) <- c("PANEL", "WGT", "ID", "CHR", "P0", "P1", "N")
      rownames(pos_file) <- c(1:nrow(pos_file))
  
      write.table(pos_file, file=paste0(wd3, 'chunk', j_num, '.pos'), col.names=T, row.names=F, quote=F, sep='\t')
      if(run_MR){
        write.table(MR_exposures, file=paste0(wd3, 'chunk', j_num, '_MRinput.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
      }
    }
    parallel::stopCluster(cl)
  }
  
#########################################
# Run TwoSample MR and Steiger filtering
#########################################

  if(run_MR){
    message("\n")
    message("\n")
    message("\n")
    message("STARTING TWO SAMPLE MR ANALYSIS...")
    message("\n")
    message("\n")
    message("\n")
    
    cl <- parallel::makeForkCluster(65, outfile = "")
    doParallel::registerDoParallel(cl)
    
    foreach (j_num = c(1:65), .combine='c', .packages = "ctwas") %dopar% {
      MR_exposures <- fread(paste0(wd3, 'chunk', j_num, '_MRinput.tsv'), stringsAsFactors = F, data.table=F)
      exposures_mr <- TwoSampleMR::format_data(MR_exposures,
                                             type = "exposure",
                                             snp_col = "SNP",
                                             beta_col = "slope",
                                             se_col = "slope_se",
                                             pval_col = "pval_nominal",
                                             log_pval=FALSE,
                                             effect_allele_col = "ALLELE1",
                                             other_allele_col = "ALLELE0",
                                             chr_col = "CHROM",
                                             phenotype_col = "ensembl_id",
                                             samplesize_col = "N")
    
      if("n_cases" %in% colnames(gwas)){
        gwas$units <- "log odds"
        outcome_mr <- TwoSampleMR::format_data(gwas[gwas$id%in%MR_exposures$SNP, ],
                                             type = "outcome",
                                             snp_col = "id",
                                             beta_col = "beta",
                                             se_col = "se",
                                             effect_allele_col = "A1",
                                             other_allele_col = "A2",
                                             eaf = "eafreq",
                                             units_col = "units",
                                             ncase_col = "n_cases",
                                             ncontrol_col = "n_controls"
        )
        outcome_mr$prevalence.outcome <- gwas[1, 'pop_prev']
      }else{
        outcome_mr <- TwoSampleMR::format_data(gwas[gwas$id%in%MR_exposures$SNP, ],
                                               type = "outcome",
                                               snp_col = "id",
                                               beta_col = "beta",
                                               se_col = "se",
                                               eaf = "eafreq", 
                                               effect_allele_col = "A1",
                                               other_allele_col = "A2",
                                               samplesize_col = "n"
        )
      }
      outcome_mr$outcome <- gwas_name
    
      dat <- TwoSampleMR::harmonise_data(
        exposure_dat = exposures_mr,
        outcome_dat = outcome_mr,
        action = 2
      )
      res <- TwoSampleMR::mr(dat)
      res <- TwoSampleMR::generate_odds_ratios(res)
      res_eggerintercept <- TwoSampleMR::mr_pleiotropy_test(dat)
      steiger <- TwoSampleMR::steiger_filtering(dat)
      
      message(paste0("Saving MR results for chunk", j_num, "..."))     
      save(res, res_eggerintercept, steiger, file=paste0(wd3, 'chunk', j_num, '_MRresults.RData'))
    }
    parallel::stopCluster(cl)
    
    res_MR_tot <- NULL
    res_eggerintercept_tot <- NULL
    res_steiger_tot <- NULL
    for(j_num in c(1:65)){
      message(paste0("Reading in MR results for chunk", j_num, "..."))
      load(paste0(wd3, 'chunk', j_num, '_MRresults.RData'))
      res_MR_tot <- rbind(res_MR_tot, res)
      res_eggerintercept_tot <- rbind(res_eggerintercept_tot, res_eggerintercept)
      res_steiger_tot <- rbind(res_steiger_tot, steiger)
      message(paste0("Deleting MR results for chunk", j_num, "..."))
      system(paste0("rm -rf ", wd3, 'chunk', j_num, '_MRresults.RData'))
      system(paste0("rm -rf ", wd3, 'chunk', j_num, '_MRinput.tsv'))
    }
    message(paste0("\n"))
    message(paste0("\n"))
    message(paste0("SAVING MR results for ALL chunks to: "))
    message(paste0(wd3, '/MRresults.RData'))
    save(res_MR_tot, res_eggerintercept_tot, res_steiger_tot, file=paste0(wd3, '/MRresults.RData'))
  }

}



