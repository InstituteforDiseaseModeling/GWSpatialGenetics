#######################################################################
# Functions to covert diploid VCF to haploid for analyses
# Author: Jessica Ribado, Institute for Disease Modeling
# Date: 2021/11 
# Update: 2024/04. This script is a duplicate of the first steps in the analysis manifest R script. Duplicating here to replicate the filtering steps to work with the NextStrain Snakemake pipeline, however the authors recommend creating the haploid VCF and associated files manually to confirm ideal filtering parameters.
#######################################################################

#######################################################################
# set-up
#######################################################################
for(p in c('vcfR', 'data.table', 'dplyr', 'tidyr')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p, repos =  "https://cloud.r-project.org", dependencies = T )
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# set global options
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE)

# load variables from snakemake
diploid_vcf <- snakemake@input[[1]]
haploid_vcf <- snakemake@output[[1]][1]
samp_missing <- snakemake@params[['samp_missing']]
site_missing <- snakemake@params[['site_missing']]
min_gt_depth <- snakemake@params[['min_read_depth']]
proportion  <- snakemake@params[['het_threshold']]

# set other relevant variables
batch_name <- gsub("_joint.*", "", basename(haploid_vcf))
output_dir <- dirname(haploid_vcf)


#######################################################################
# functions
#######################################################################
# VCF files
filter_batch_duplicates <- function(samples){
  # return duplicates of files between batches for filtering 
  duplicated_samples <- c(gsub(".batch.*", "",  
    samples[grepl("\\.batch", samples)]),
    samples[grepl("MISC|BAT", samples)])
  return(duplicated_samples)  
}


vcf2haploid <- function(vcf){ 
  #  vcf <- vcfR::read.vcfR(vcf_files[[1]], verbose = FALSE)
  tmp_vcf <- vcf
  # reset heterozygous genotypes as a single representation
  gt <- extract.gt(tmp_vcf)
  # tmp_vcf@gt[,-1][is_het(gt)] <- "N"
  is.na(tmp_vcf@gt[,-1][is_het(gt)]) <- TRUE
  het_gt <- extract.gt(tmp_vcf)
  het_gt <- data.frame(apply(het_gt, 2, function(x) gsub("\\|", "\\/", x)))
  het_hap <- data.frame(apply(het_gt, 2, function(x) gsub("\\/.*", "", x)))

  # revert heterozygous calls to homozygous if majority of calls go to one allele
  tmp_vcf <- vcf
  gt <- extract.gt(tmp_vcf)
  is.na(tmp_vcf@gt[,-1][!is_het(gt)]) <- TRUE
  # Extract allele depths.
  ad <- extract.gt(tmp_vcf, element = "AD")
  ad1 <- masplit(ad, record = 1)
  ad2 <- masplit(ad, record = 2)
  freq1 <- ad1/(ad1+ad2)
  # create a copy of the matrix to convert to homozygous
  het_hap_adj <- het_hap
  het_hap_adj[freq1 >= proportion] <- "0"
  het_hap_adj[freq1 <= 1-proportion] <- "1"
  
  # remove genotypes with less than the minimum read depth
  tmp_vcf <- vcf
  total_dp <- extract.gt(tmp_vcf, element = 'DP', as.numeric = TRUE)
  is.na(het_hap_adj) <- is.na(total_dp) <- total_dp < min_gt_depth
  is.na(het_hap_adj[het_hap_adj == "\\*"]) <- TRUE

  # paste with other sample information
  gt2 <- extract.gt(tmp_vcf, extract = FALSE)
  is.na(gt2) <- is.na(het_hap_adj)
  gt_merge <- matrix(paste(as.matrix(het_hap_adj), gt2, sep=":"), nrow=nrow(gt), dimnames=dimnames(gt))
  is.na(gt_merge[gt_merge == "NA:NA"]) <- TRUE
  tmp_vcf@gt[,-1] <- gt_merge
  
  # site filters
  # filter out sites with a majority proportion of missing reads
  gt_counts <- cbind.data.frame(
    ref =  rowSums(het_hap_adj == "0", na.rm=T),
    missing = rowSums(is.na(het_hap_adj))
  )
  gt_counts$alt = ncol(het_hap_adj) - rowSums(gt_counts)
  
  gt_counts <- dplyr::filter(gt_counts, alt > 0 & ref > 0) 
  
  site_keep_sing <- dplyr::filter(gt_counts, missing < site_missing * nrow(gt_counts))
  site_keep <- dplyr::filter(gt_counts, alt > 1 & ref > 1 & missing < site_missing * nrow(gt_counts))
  
  het_site_rm <- het_hap_adj[row.names(het_hap_adj) %in% row.names(site_keep_sing), ]
  
  # sample filters
  # identify duplicate samples 
  dup_samples <- filter_batch_duplicates(names(het_hap_adj))
  # check for duplicates that show up more than once
  dup_removal <- names(table(dup_samples))[table(dup_samples) <= 1]
  double_dups <- names(table(dup_samples))[table(dup_samples) > 1]
  if (length(double_dups) > 0){
    for (duplicate in double_dups){
      duplicates <- names(het_hap_adj)[grepl(duplicate, names(het_hap_adj))]
      dup_counts <- colSums(is.na(het_hap_adj[,duplicates]))
      dup_removal <- c(dup_removal, names(dup_counts)[dup_counts < max(dup_counts)])
    }
  }
  duplicated_keep <- !names(het_site_rm) %in% dup_removal

  # identify samples with a high proportion of positions not meeting minimum read depth  
  missing_percent <- colSums(is.na(het_site_rm))/nrow(het_site_rm)
  depth_keep <- missing_percent < samp_missing
  # save read depths to report back to wetlab for potential edge cases 
  missing_df <- data.table(lab_id=names(missing_percent[duplicated_keep]), filt_site_missing_prop=missing_percent[duplicated_keep]) %>%
    dplyr::mutate(excluded_for_analysis=ifelse(filt_site_missing_prop < samp_missing, "No", "Yes"))
  
  samp_keep <- duplicated_keep & depth_keep

  # update the VCF object  
  # remove batch information for duplicated samples
  tmp_vcf@gt <- tmp_vcf@gt[,c(TRUE, samp_keep)]
  # remove batch information for duplicated samples
  colnames(tmp_vcf@gt) <- gsub("-batch.*" , "", colnames(tmp_vcf@gt))
  colnames(tmp_vcf@gt) <- gsub("-" , ".", colnames(tmp_vcf@gt))
  tmp_vcf <- tmp_vcf[row.names(gt) %in% row.names(site_keep), ]

  # organize summary statistics for filtering
  summary_df <- data.frame(rbind(
    "Diploid VCF  " = basename(diploid_vcf),
    "Samples (All) " = ncol(vcf@gt),
    "Samples (Exluding technical replicates) " = ncol(vcf@gt) - length(dup_samples),
    "Variants " = nrow(het_hap),
    "Minimum genotype ready depth " = min_gt_depth,
    "Heterozygous recode proportion " = paste(1 - proportion, proportion, sep=","),
    "Maximum proportion of samples missing a call " = site_missing,
    "Maximum proportion of missing sites per sample " = samp_missing,
    "Samples after filtering  " = sum(samp_keep),
    "Variants after filtering " = nrow(site_keep_sing),
    "Variants after filtering, singletons excluded  " = nrow(site_keep)
  )) %>% tibble::rownames_to_column()
  names(summary_df) <- c("summary", "value")
  
  return(list(summary = summary_df,
              vcf = tmp_vcf,
              missing_df = missing_df))
}  

# other functions
SaveTabDelim <- function(r_obj, path){
  write.table(r_obj, file = path, quote = F, row.names = F, sep="\t")
}


################################################################################
# format data files
################################################################################
# convert diploid VCF to haploid, filter variants and samples that are lower quality
vcf <- vcfR::read.vcfR(diploid_vcf, verbose = FALSE)
hap_vcf <- vcf2haploid(vcf)

# save files for other analyses - i.e. Nextstrain, for collaborators
vcfR::write.vcf(hap_vcf$vcf,  paste0(output_dir, "/", batch_name, "_jointHaploidFilter.vcf.gz"))
SaveTabDelim(hap_vcf$summary, paste0(output_dir, "/", batch_name, "_filterSummary.tsv"))
