#######################################################################
# Functions to covert diploid VCF to haploid for analyses
# Author: Jessica Ribado, Institute for Disease Modeling
# Date: 2021/11 
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

#######################################################################
# functions
#######################################################################
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
  
  # apply filter
  # filter samples with low read depth 
  samp_keep <- (colSums(is.na(het_hap_adj))/nrow(het_hap_adj) < samp_missing) & 
    sapply(names(het_hap_adj), function(y) !grepl(y, "BAT"))
  het_samp_rm <- het_hap_adj[, samp_keep]
  
  # filter out sites with a majority proportion of missing reads
  gt_counts <- cbind.data.frame(
     ref =  rowSums(het_samp_rm == "0", na.rm=T),
     missing = rowSums(is.na(het_samp_rm))
   )
  gt_counts$alt = ncol(het_samp_rm) - rowSums(gt_counts)
   
  site_keep_sing <- dplyr::filter(gt_counts, alt > 0 & ref > 0 & missing < site_missing * ncol(het_hap_adj))
  site_keep <- dplyr::filter(gt_counts, alt > 1 & ref > 1 & missing < site_missing * ncol(het_hap_adj))
  
  tmp_vcf@gt <- tmp_vcf@gt[, c(TRUE, samp_keep)]
  tmp_vcf <- tmp_vcf[rownames(het_hap_adj) %in% row.names(site_keep), ]

  
  # organize summary statistics for filtering 
  summary_df <- data.frame(rbind(
    "Diploid VCF" = basename(diploid_vcf),
    "Samples" = ncol(het_hap),
    "Variants" = nrow(het_hap),
    "Minimun genotype ready depth" = min_gt_depth,
    "Heterozygous recode proportion" = paste(1 - proportion, proportion, sep=","),
    "Maximum proportion of samples missing a call" = site_missing,
    "Maximum proportion of missing sites per sample" = samp_missing,
    "Haploid VCF" = basename(haploid_vcf),
    "Samples after filtering" = sum(samp_keep),
    "Variants after filtering" = nrow(site_keep_sing),
    "Variants after filtering, singletons excluded" = nrow(site_keep)
  )) %>% tibble::rownames_to_column() 
  names(summary_df) <- c("summary", "value")
  
  return(list(summary = summary_df,
              vcf = tmp_vcf))
}  

#######################################################################
# run and save output
#######################################################################
vcf <- vcfR::read.vcfR(diploid_vcf, verbose = FALSE)
filt_vcf <- vcf2haploid(vcf)
write.table(filt_vcf$summary, file = gsub("_jointHaploidFilter.vcf.gz", "_filterSummary.tsv", haploid_vcf), 
sep="\t", row.names = F, col.names=F, quote=F)
vcfR::write.vcf(filt_vcf$vcf, haploid_vcf)
