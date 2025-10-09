################################################################################
# mtDNA Recalibration with complete Set Quick Check 
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: April 2024
# Updated August 2025
################################################################################


################################################################################
# Characterize BQSR contributing samples - Manual run 
################################################################################
# ids <- data.table::fread("/mnt/data/guinea_worm/bsqr_hq_samples.txt", header = F)
# 
# # Regex to extract parts: 3 letters, 3 letters, 4 digits (year)
# matches <- str_match(ids$V1, "^([A-Z]{3})([A-Z]{3})?(\\d{4})")
# 
# # Create data frame
# result <- data.frame(
#   full_id = ids,
#   sample_id = gsub("-.*", "", ids$V1),
#   country_code = matches[,2],
#   host_code = matches[,3],
#   year = matches[,4],
#   stringsAsFactors = FALSE
# ) %>%
#   dplyr::filter(country_code != "BAT")
# unique(result$sample_id) %>% length()    # N specimens
# table(result$country_code, result$year)  # Country specimen counts


################################################################################
# Load data
################################################################################
original <- readRDS("/mnt/data/guinea_worm/analyses/batch_Apr042025/GenomicSamplesMetadata_Database_batch_Apr042025.rds")
new <- readRDS("/mnt/data/guinea_worm/analyses/recalJuly2025/GenomicSamplesMetadata_Database_recalJuly2025.rds")

################################################################################
# Overlap of barcodes
################################################################################
positions_old <- read.delim("/mnt/data/guinea_worm/analyses/batch_Apr042025/batch_Apr042025_BarcodeVariantSummary.tsv", sep="\t")
positions_new <- read.delim("/mnt/data/guinea_worm/analyses/recalJuly2025/recalJuly2025_BarcodeVariantSummary.tsv", sep="\t")

overlap <- intersect(positions_old$POS, positions_new$POS)
print(paste("Overlapping barcode variant positions", length(overlap)))
print(paste("Original unique barcode variants positions", nrow(positions_old) - length(overlap)))
print(paste("New unique barcode variant positions", nrow(positions_new) - length(overlap)))


################################################################################
# Unique amplicon counts
################################################################################
original <- readRDS("/mnt/data/guinea_worm/analyses/batch_Apr042025/GenomicSamplesMetadata_Database_batch_Apr042025.rds")
new <- readRDS("/mnt/data/guinea_worm/analyses/bsqr2/GenomicSamplesMetadata_Database_bsqr2.rds")

amplicon_counts <- inner_join(
  original %>% group_by(country) %>% summarize(original = n_distinct(amplicon_barcode)),
  new %>% group_by(country) %>% summarize(new = n_distinct(amplicon_barcode))
)


################################################################################
# Check siblings
################################################################################
siblings <- read.delim("/mnt/data/guinea_worm/analyses/dcifer_test/trial_dataset/SSU_SibPairingsForJessica_23.12.06microsatdata.txt", sep="\t") %>%
  dplyr::filter(`FSBinary...0.97.` == 1)
sibling_names <- unique(c(siblings$PutativeSibID1, siblings$PutativeSibID2))

siblings_old_barcodes <- dplyr::filter(original, genomics_sample_id %in% sibling_names) %>% dplyr::select(genomics_sample_id, amplicon_barcode)
siblings_new_barcodes <- dplyr::filter(new, genomics_sample_id %in% sibling_names) %>% dplyr::select(genomics_sample_id, amplicon_barcode)
merged_bc <- 
  dplyr::left_join(siblings, dplyr::rename(siblings_new_barcodes, 'PutativeSibID1'='genomics_sample_id', 'new1'= 'amplicon_barcode')) %>%
  dplyr::left_join(., dplyr::rename(siblings_new_barcodes, 'PutativeSibID2'='genomics_sample_id', 'new2'= 'amplicon_barcode')) %>%
  dplyr::left_join(., dplyr::rename(siblings_old_barcodes, 'PutativeSibID1'='genomics_sample_id', 'old1'= 'amplicon_barcode')) %>%
  dplyr::left_join(., dplyr::rename(siblings_old_barcodes, 'PutativeSibID2'='genomics_sample_id', 'old2'= 'amplicon_barcode')) %>%
  na.exclude %>%
  dplyr::mutate(same_new = ifelse(new1 == new2, 1, 0),
                same_old = ifelse(old1 == old2, 1, 0))

# check pairwise
t <- dplyr::filter(barcode_pairwise, amplicon_barcode_conservative_pair %in% c("781_782", "781_783", "782_783"))




################################################################################
# Check Chad clusters
################################################################################
amdabri_specimens <- c('PDB18-047', 'PDB18-093_L1', 'PDB18-093', 'PDB19-001', 'PDB19-124', 'PDB21-016')
amhabile_specimens <- c("PDB18-048", "PDB18-055", "PDB18-056", "PDB18-072", "PDB18-073", "PDB19-043", "PDB19-046", "PDB19-125") #, "12")

require('optparse')
option_list = list(
  make_option(c("-f", "--filt_prop"), type="numeric", default=0.1, 
              help="The proportion of missing positions, after filtering high quality variants, tolerated in a barcode for analyses. Default is 10% of the barcode length.", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

output_dir <- '/mnt/data/guinea_worm/analyses/recalJuly2025/'
script_dir <- ('/home/jribado/git/GWSpatialGenetics/analyses')
source(paste(script_dir, "mtDNA_common_functions.R", sep="/"))
source(paste(script_dir, "mtDNA_plotting.R", sep="/"))
barcode_pairwise <- data.table::fread('/mnt/data/guinea_worm/analyses/recalJuly2025/recalJuly2025_relatednessBarcodes.txt')
#barcode_pairwise <- data.table::fread('/mnt/data/guinea_worm/analyses/batch_Apr042025/batch_Apr042025_relatednessBarcodes.txt')

ClusterPariwiseMatrices <- function(cluster_name, df = metadata){
  cluster_df <- dplyr::filter(df, cluster == !!cluster_name & analysis_inclusion == "Included") 
  if(nrow(cluster_df) > 2){
    cluster_diff <- SubsampleMerge(cluster_df, split_variable = NULL) %>%
      mutate(missing = missing/n_variants)
    ClusterPlots(df = cluster_df, pairwise_df = cluster_diff,
                 #sample_name="wormnum_dpdxid",
                 output_name=paste0("cluster_", gsub(" |cluster", "", cluster_name)))
  }
}

small_df <- dplyr::filter(new, wormnum_dpdxid %in% c(amdabri_specimens, amhabile_specimens) | genomics_sample_id == 'CHDHUM2019_05377') %>%
  dplyr::mutate(cluster = ifelse(wormnum_dpdxid %in% amdabri_specimens, "Amdabri", "Amhabile"))
missing_count_variable = names(new)[grepl("missing_n_", names(new))]
n_variants <- readr::parse_number(missing_count_variable)

for (j in c("Amdabri", "Amhabile")){
  ClusterPariwiseMatrices(j, small_df)
}

