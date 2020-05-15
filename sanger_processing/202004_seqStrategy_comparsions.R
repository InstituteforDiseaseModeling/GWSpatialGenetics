# Rarefactions for Guinea worm sequencing 
# Author: Jessica Ribado
# Date: 01/24/2020

###############################################################################
# setup
############################################################################### 
# load libraries
for(p in c('data.table', 'dplyr', 'tidyverse',
           'protr', 'Biostrings', 'doParallel',
           'geosphere', 'lubridate',
           'ggplot2', 'ggforce')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE)


source("~/git/DDA-Genetics-GuineaWormTransmissionInChad/ngs_processing/202004_seqStrategy_functions.R")
theme_set(theme_j())

############################################################################### 
# specify data directories and files
###############################################################################
# project directory
gw_parent <- "/home/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Guinea Worm Genetics in Chad/gw_bwaAlign"
gw_dir <- paste(gw_parent, "ngs_analyses", sep="/")

# data files
mt_gatk_file <- "/mnt/md0/guinea_worm/ngs_analyses/unmatched/mt/01_processing/03_variant_calls/jointGenotype_2Iter_filtered.txt"
mt_bcf_dir   <- '/mnt/md0/guinea_worm/ngs_analyses/unmatched/mt/01_processing/03_variant_calls/bcftools/individual_vcfs'
mt_bcf_files <- list.files(mt_bcf_dir, pattern=".vcf.gz$")


################################################################################
# whole mitocondrial distribution compared to mt gene distribution for similarity
################################################################################ 
# mnaually remove regions with heterozygous calls
het_sites <- c(seq(5950, 6670)) 
# mt_gatk <- read.delim(mt_gatk_file)
mt_gatk <- var_group(mt_gatk_file) %>% dplyr::filter(!POS %in% het_sites)
mt_gatk <- cbind(mt_gatk[,c(1:4)], gt_count(mt_gatk[,-c(1:4)]))

# calculate similairty
gatk_missingRm <-dplyr::filter(mt_gatk, unique_bases > 1, N < 2, DP_min == 0)#, singleton_bases > 1) 
gatk_sim <- rbind.data.frame(
  sim_df(gatk_missingRm) %>% dplyr::mutate(region="All loci"),
  sim_df(dplyr::filter(gatk_missingRm, bc_region == "True")) %>% dplyr::mutate(region="Within loci"), 
  sim_df(dplyr::filter(gatk_missingRm, bc_region == "False")) %>% dplyr::mutate(region="Outside loci")
) 

gatk_countT <- table(gatk_missingRm$bc_region)[['True']]
gatk_countF <- table(gatk_missingRm$bc_region)[['False']]
pairwise_plots(gatk_sim, gatk_countT, gatk_countF)
# run permutation tests for the one above
set.seed(123)
reps <- 50
sim_boot <- do.call("rbind", lapply(rep(1:reps), function(i){
  print(i)
  false_sub <- dplyr::filter(gatk_missingRm, bc_region =="False") %>% 
    dplyr::sample_n(gatk_countT)
  mt_df <- sim_df(false_sub) 
  mt_df$region = paste0("Outside loci subset\n(n=", gatk_countT, ")")
  mt_df$boot <- i
  return(mt_df)
}))  
# merge in the full set to combine in one figure
gatk_simScores <- rbind(dplyr::mutate(gatk_sim, boot = 0), sim_boot,
                        dplyr::filter(gatk_sim, region=="Within loci") %>% 
                          dplyr::mutate(boot = 51)) %>%
  dplyr::mutate(facet = ifelse(boot == 0, "All", "Subset"))


# plot 
sim_pdf <- gatk_simScores %>% 
  ggplot(aes(x=value, fill=region, group=paste(region, boot))) + 
  geom_histogram(aes(y = ..density..), bins=20, alpha=0.2, position="identity") +
  geom_density(alpha=0.5) +
  scale_fill_manual(values = c("pink", "black", "grey80", "red"),
                    labels = c(paste0("All variants\n(n=", gatk_countT + gatk_countF, ")"),
                               paste0("Outside loci\n(n=", gatk_countF, ")"),
                               paste0("Outside loci subsample\n(n=", gatk_countT, ")"),
                               paste0("Within loci\n(n=", gatk_countT, ")")),
                    name = "") + 
  labs(x="Pairwise genetic similarity", y="Comparisons") + 
  theme(legend.position = "top") +
  facet_grid(~facet) +
  xlim(c(0,1))

ks.test(dplyr::filter(gatk_sim, region == "Within loci") %>% .[["value"]],
        dplyr::filter(gatk_sim, region == "Outside loci") %>% .[["value"]])
ks.test(dplyr::filter(gatk_sim, region == "Within loci") %>% .[["value"]],
        dplyr::filter(gatk_sim, region == "All loci") %>% .[["value"]])
ggsave("20200502_seqStrategy_gatkParams_rm5950-6670.png", plot = sim_pdf, path = gw_dir, width = 5.5, height = 4, units = c("in"), dpi = 300)
#ggsave("20200502_seqStrategy_gatkParams_rm5950-6670_singEx.png", plot = sim_pdf, path = gw_dir, width = 5.5, height = 4, units = c("in"), dpi = 300)



################################################################################
# GATK v bcftool mpileup genotype calling
################################################################################ 
# read in vcf_files
vcf_List <- lapply(mt_bcf_files, function(i){
  print(i)
  mpile_vcf <- vcfR::read.vcfR(paste(mt_bcf_dir, i, sep="/"))
  vcf_df    <- cbind(as.data.frame(vcfR::getFIX(mpile_vcf)), vcfR::INFO2df(mpile_vcf), gt = mpile_vcf@gt[,2]) %>%
    dplyr::filter(is.na(IDV)) %>%
    dplyr::mutate(GT = ifelse(gt == "0", REF, ifelse(grepl("1", gt), ALT, "."))) %>%
    dplyr::select(CHROM, POS, REF, DP, GT)
  dp <- paste0(gsub("\\.vcf.gz", "", i), "_DP")
  gt <- paste0(gsub("\\.vcf.gz", "", i), "_GT")
  names(vcf_df) <- c("CHROM", "POS", "REF", dp, gt)
  return(vcf_df)
})
vcf_merge_df <- vcf_List %>% purrr::reduce(left_join, by = c("CHROM", "POS", "REF")) 

# identify positions with a variant
var_pos <- vcf_merge_df %>%
  dplyr::mutate(bc_region = ifelse(POS %in% seq(3788, 4534) | POS %in% seq(2628, 3345) | POS %in% seq(12562, 14566), "True", "False"),
                POS = as.numeric(POS)) %>%
  dplyr::filter_at(vars(ends_with("_DP")), all_vars(. > 10)) %>%
  dplyr::select(CHROM, POS, bc_region, ends_with("_GT")) 

# count the number of bases at each position
var_comb <- cbind(var_pos[,c(1:3)], gt_count(var_pos[,-c(1:3)]))

################################################################################ 
# functions to run similarity scores
sim_scores <- function(pos_df){
  df_t <- tidyr::unite(data.frame(t(pos_df[,c(4:20)])), barcode, seq(1, nrow(pos_df)), sep="")
  var_barcodes <- as.list(df_t[,1])
  names(var_barcodes) <- row.names(df_t)
  var_mat <- simMat(var_barcodes, 8, gsub("_GT", "", names(var_barcodes)))
  var_df  <- simMat2df(var_mat)
  var_df$bc_length <- nchar(var_barcodes[[1]])
  return(var_df)
}

sim_all <- rbind.data.frame(
  sim_scores(dplyr::filter(var_comb, bc_region=="True", unique_bases != 1)) %>% dplyr::mutate(region="True"),
  sim_scores(dplyr::filter(var_comb, bc_region=="False", unique_bases != 1)) %>% dplyr::mutate(region="False"))
sim_singEx <- rbind.data.frame(
  sim_scores(dplyr::filter(var_comb, bc_region=="True", singleton_bases > 1)) %>% dplyr::mutate(region="True"),
  sim_scores(dplyr::filter(var_comb, bc_region=="False", singleton_bases > 1)) %>% dplyr::mutate(region="False"))

################################################################################ 
# plotting
pairwise_plots <- function(sim_df){
  ggplot(sim_df, aes(x=value, group=region, fill=region)) + 
    geom_histogram(aes(y = ..density..), bins=20, alpha=0.2, position="identity") +
    geom_density(alpha=0.5) +
    scale_fill_manual(values = c("black", "red"),
                       labels = c(paste0("Outside of loci\n(n=", unique(sim_df$bc_length[which(sim_df$region == 'False')]), ")"),
                                  paste0("Within loci\n(n=", unique(sim_df$bc_length[which(sim_df$region == 'True')]), ")")),
                       name = "") + 
    labs(x="Pairwise genetic similarity", y="Comparisons") + 
    theme(legend.position = "top") +
    xlim(c(0,1))
}

ggsave("202005_seqStrategy_bcftoolsVariantCheck.png", plot = pairwise_plots(sim_all), path = gw_dir, width = 5.5, height = 4, units = c("in"), dpi = 300)  
ggsave("202005_seqStrategy_bcftoolsVariantCheck_singEx.png", plot = pairwise_plots(sim_singEx), path = gw_dir, width = 5.5, height = 4, units = c("in"), dpi = 300)   


################################################################################
# GATK v bcftool mpileup genotype positions comparison
################################################################################  
# combine coutns from both technologies
gatk_variants <- dplyr::filter(gatk_missingRm, unique_bases > 1) %>%
  dplyr::select(CHROM, POS, REF, ALT, A, C, G, T) 
bcftools_variants <- dplyr::filter(var_comb, unique_bases > 1) %>%
  dplyr::select(CHROM, POS, A, C, G, T) 
pipeline_variants <- dplyr::full_join(gatk_variants, bcftools_variants, by=c("CHROM", "POS"))
pipeline_variants$shared_genotypes <- unlist(lapply(seq(1, nrow(pipeline_variants)), function(j){
  all(pipeline_variants[j,5:8] == pipeline_variants[j,9:12])
}))  
write.table(pipeline_variants, paste(gw_dir, "variantCaller_counts.txt", sep="/"), sep="\t", quote =  F, row.names = F)

pipeline_positions <- dplyr::filter(pipeline_variants, shared_genotypes == "TRUE") %>% .[["POS"]]
gatk_shared <-dplyr::filter(mt_gatk, POS %in% pipeline_positions)
gatk_simShared <- rbind.data.frame(
  sim_df(gatk_shared) %>% dplyr::mutate(region="All loci"),
  sim_df(dplyr::filter(gatk_shared, bc_region == "True")) %>% dplyr::mutate(region="Within loci"), 
  sim_df(dplyr::filter(gatk_shared, bc_region == "False")) %>% dplyr::mutate(region="Outside loci")
) 

gatk_countT <- table(gatk_shared$bc_region)[['True']]
gatk_countF <- table(gatk_shared$bc_region)[['False']]
shared_dist <- pairwise_plots(gatk_simShared, gatk_countT, gatk_countF)
ggsave("202005_seqStrategy_sharedVar.png", plot = shared_dist, path = gw_dir, width = 5.5, height = 4, units = c("in"), dpi = 300)  
