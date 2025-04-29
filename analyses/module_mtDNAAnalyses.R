#!/usr/bin/env Rscript
################################################################################
# mtDNA Processing and cross-country analysis pipeline 
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: April 2024
# Updated February 2025
################################################################################

################################################################################
# input files and options
################################################################################
# set command line and default options
require('optparse')
option_list = list(
  make_option(c("-d", "--diploid_vcf"), type="character", default=NULL, 
              help="Diploid VCF generated from the Guinea worm NGS tiling amplicon panel processing pipeline.", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Metadata for the sequenced samples in the VCF file. Minimum fields include sampleID (Vassar generated), country, year, and host for standard analyses.", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=paste0(format(Sys.time(), "%Y%m%d"), "_run"), 
              help="Output directory. A subdirectory with the batch name from the VCF file name will be automatically created at the provided path. If none provided, will output to a folder in the running directory with the date.", metavar="character"),
  make_option(c("-s", "--samp_missing"), type="numeric", default=0.5, 
              help="The proportion of positions missing in a sample to consider sample exclusion.", metavar="character"),
  make_option(c("-p", "--site_missing"), type="numeric", default=0.7, 
              help="The proportion of samples missing a variant call to consider variant exclusion from the barcode.", metavar="character"),
  make_option(c("-g", "--min_gt_depth"), type="numeric", default=5, 
              help="The minimum number of reads per site in individual sample to consider the variant call accurate.", metavar="character"),
  make_option(c("-e", "--het_proportion"), type="numeric", default=1, 
              help="The proportion cut off for heterozygous calls (i.e. from heteroplasmy or potential nuclear genome insertions). Default parameter is 1 which removes any heterozygous calls and report that position as missing.", metavar="character"),
  make_option(c("-f", "--filt_prop"), type="numeric", default=0.1, 
              help="The proportion of missing positions, after filtering high quality variants, tolerated in a barcode for analyses. Default is 10% of the barcode length.", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


################################################################################
# set-up
################################################################################
# load libraries
packages_to_install <- c('Rcpp', 'parallel', 'vcfR', 'geodist', 'lubridate',
                         'data.table', 'dplyr', 'tidyr', 'tibble', 'stringr', 
                         'ggplot2', 'ggpubr', 'iNEXT')
for(p in packages_to_install){
  if(!p %in% installed.packages()[,1]){
    install.packages(p, repos = "http://cran.us.r-project.org") 
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}  

# load functions
# function to get the running script path if running in another directory
# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

# script_dir <- getScriptPath()
setwd('/home/jribado/git/GWSpatialGenetics/analyses')
script_dir <- try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
source(paste(script_dir, "module_checkMetadata.R", sep="/"))
source(paste(script_dir, "module_vcfdiploid2haploid.R", sep="/"))
source(paste(script_dir, "module_barcodeFunctions.R", sep="/"))
source(paste(script_dir, "mtDNA_common_functions.R", sep="/"))
source(paste(script_dir, "mtDNA_plotting.R", sep="/"))


# plotitng options
theme_set(theme_bw())

diploid_vcf = '/mnt/data/guinea_worm/processing/vcf_files/batch_Apr042025_jointGenotypeFiltered.vcf.gz'
metadata_file = '/mnt/data/guinea_worm/metadata/GenomicSamplesMetadata_Database_v4.2_250326.rds'

# update parser if running in manual mode
if(is.null(opt$diploid_vcf)){
  opt$diploid_vcf <- diploid_vcf
  opt$metadata <- metadata_file
  opt$output_dir <- "/mnt/data/guinea_worm/analyses"
}

batch_name <- gsub("_jointGenotype.*", "", basename(opt$diploid_vcf))
output_dir <- paste(opt$output_dir, batch_name, sep="/")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

# confirm VCF and metadata files exist 
if(!file.exists(opt$diploid_vcf) | is.null(opt$diploid_vcf)){
  stop("VCF file not provided or user defined file does not exist. Confirm the file path and try again.")
}

if(!file.exists(opt$metadata) | is.null(opt$metadata)){
  stop("Metadata file not provided or user defined file does not exist. Confirm the file path and try again.")
}


################################################################################
# format data files
################################################################################
# convert diploid VCF to haploid, filter variants and samples that are lower quality
vcf <- vcfR::read.vcfR(opt$diploid_vcf, verbose = FALSE)
hap_vcf <- vcf2haploid(vcf)

# generate barcodes from high quality positions
vcf_clust <- gt2barcode(hap_vcf$vcf) 

# get names of samples in the vcf
vcf_samples <- colnames(hap_vcf$vcf@gt)[-1]
vcf_numbers <- cbind.data.frame(
  sample = vcf_samples,
  vassar_worm = gsub(".*_", "", vcf_samples))

# read in and check metadata
if(grepl("\\.txt$", opt$metadata)){
  metadata <- data.table::fread(opt$metadata, sep="\t")
} else if(grepl("\\.rds$", opt$metadata)){
  metadata <- readRDS(opt$metadata)
} else{
  stop("Metadata does not have .txt or .rds extension. User must manually code data input.")
}
metadata[metadata == "." & !is.na(metadata)] <- NA
metadata <- run_metadata_checks(metadata, 
                                vcf_samples = vcf_numbers, 
                                excluded_samples = hap_vcf$missing_df)
metadata <- dplyr::left_join(metadata, vcf_clust) %>%
  dplyr::mutate(epi_foci = ifelse(epi_foci == "", NA, gsub(" focal area", "", epi_foci)))


# save VCF files and filter for other analyses 
vcfR::write.vcf(hap_vcf$vcf,  paste0(output_dir, "/", batch_name, "_jointHaploidFilterWithSingletons.vcf.gz"))
SaveTabDelim(hap_vcf$summary, paste0(output_dir, "/", batch_name, "_filterSummary.tsv"))
SaveTabDelim(data.frame(hap_vcf$vcf@fix), paste0(output_dir, "/", batch_name, "_BarcodeVariantSummary.tsv"))
rm(vcf)


################################################################################
# Genetic similarity
################################################################################
# Only compare the unique barcodes to reduce compute and intermediate file sizes
unique_sequences <- names(sort(table(metadata[['sequence']]), decreasing = T))
barcode_pairwise <- BarcodePairwiseParallel(unique_sequences)
barcode_pairwise <- barcode_pairwise %>%  
  dplyr::inner_join(., unique(dplyr::select(metadata, sequence, amplicon_barcode_conservative)) %>% dplyr::rename(sequence.x=sequence), by="sequence.x") %>%
  dplyr::inner_join(., unique(dplyr::select(metadata, sequence, amplicon_barcode_conservative)) %>% dplyr::rename(sequence.y=sequence), by="sequence.y") 

barcode_pairwise <- dplyr::rowwise(barcode_pairwise) %>%
  dplyr::mutate(amplicon_barcode_conservative_pair = 
                  paste(min(amplicon_barcode_conservative.x, amplicon_barcode_conservative.y), 
                        max(amplicon_barcode_conservative.x, amplicon_barcode_conservative.y), sep="_"))
SaveTabDelim(barcode_pairwise, paste0(output_dir, "/", batch_name, "_relatednessBarcodes.txt"))


################################################################################
# Manually cluster barcodes with few missing positions to complete barcodes
################################################################################
missing_count_variable = names(metadata)[grepl("missing_n_", names(metadata))]
n_variants <- readr::parse_number(missing_count_variable)

barcodes <- unique(dplyr::select(vcf_clust, amplicon_barcode_conservative, sequence, frequency_conservative, eval(missing_count_variable))) %>%
  dplyr::filter(!is.na(sequence)) 
complete_barcodes <- dplyr::filter(barcodes, get(missing_count_variable) == 0) %>% .[['amplicon_barcode']]
updated_barcodes <- ManualBarcodeRecode(barcodes) 
# update metadata to include new non-ambiguous groupings
metadata <- dplyr::inner_join(metadata, updated_barcodes) %>%
  dplyr::mutate(complete_regroup = ifelse(amplicon_barcode %in% complete_barcodes, "Yes", "No")) %>%
  AssignAmpliconGrouping(.) %>% ungroup()

# update color for amplicons 
common_amplicons <- readr::parse_number(
  unique(metadata$amplicon)[!unique(metadata$amplicon) %in% names(reduced_base)]) 
common_amplicons <- unique(common_amplicons[!is.na(common_amplicons) & common_amplicons>0])
if(length(common_amplicons) <= length(all_colors)){
  print("Sufficient colors in default pallete for plotting barcode lineages.")
  all_colors <- setNames(base_colors2[1:length(common_amplicons)], sort(common_amplicons))
  nextstrain_colors <- c(all_colors, reduced_base)  
} else{
  print("Too many barcodes provided for default color scheme. 
        Some barcodes will not be represented in the barcode legend.")
}
metadata$amplicon <- factor(metadata$amplicon, levels=names(nextstrain_colors))
metadata$host <- factor(metadata$host, levels=names(host_colors))
metadata$amplicon_barcode <- as.numeric(metadata$amplicon_barcode)
metadata <- as.data.table(metadata)
metadata <- metadata[, gps_id := .GRP, by = .(gps_e, gps_n)]

# save metadata file with barcode information 
rm(barcodes)
SaveTabDelim(metadata, paste0(output_dir, "/", "GenomicSamplesMetadata_Database_", batch_name, ".tsv"))
saveRDS(metadata, paste0(output_dir, "/", "GenomicSamplesMetadata_Database_", batch_name, ".rds"))


################################################################################
# Visualize variant and barcode quality 
################################################################################
# barcode changes post manual grouping
bc_changes <- BarcodeClusterDifferences(metadata)

# major allele frequency distribution
samples_remove <- dplyr::filter(vcf_clust, analysis_inclusion == "Excluded") %>% .[["sample"]]
SavePlots(AlleleFrequencyPlot(hap_vcf$vcf, samples_remove), 
          output_dir, paste0(batch_name, "_MajorAlleleFrequencies.png"), 
          height=4, width=7)

# missing positions to barcode frequency relationship 
SavePlots(MissingFrequencyPlot(bc_changes), 
          output_dir, paste0(batch_name, "_ConservativeMissingFrequency.png"), 
          height=4, width=6)

# visualize confidence of barcode groups and changes post clustering
SavePlots(BarcodeChangeMulti(bc_changes), 
          output_dir, paste0(batch_name, "_BarcodeClusteringChanges.png"), 
          height=5, width=10)

rm(bc_changes, hap_vcf, vcf_clust)


################################################################################
# Barcode summaries
################################################################################
# overall barcode distribution
SavePlots(BarcodeByCountryPlot(metadata), 
          output_dir, paste0(batch_name, "BarcodeHighCaseCountries.png"), 
          height=4, width=10)

# by country barcodes
barcode_countries <- names(table(metadata$country))[table(metadata$country) > 5]
barcode_countries <- barcode_countries[barcode_countries != ""]
lapply(barcode_countries, function(i){
  BarcodeCountPlots(dplyr::filter(metadata, country == !!i))
  BarcodeCountPlots(dplyr::filter(metadata, country == !!i), column="month")
})

# barcodes by country overlap - groups mitochondrial lineages across long ranges 
count_barcodes <- dplyr::group_by(metadata, country, amplicon_barcode) %>%
  dplyr::summarise(country_count = n())
count_barcodes_total <- table(count_barcodes$amplicon_barcode)
shared_barcodes <- names(count_barcodes_total)[count_barcodes_total > 1]
potential_country_links <- unique(dplyr::bind_rows(
  lapply(shared_barcodes, Barcode2PotentialKinship)
  )) 
SaveTabDelim(potential_country_links, paste0(output_dir, "/", batch_name, "_CountryBarcodeOverlapCounts.tsv"))

# identify the specimens with shared mitochondrial DNA across countries for sibling prediction
potential_links_specimens <- dplyr::select(potential_country_links, amplicon_barcode, country1, country2, overlapping_year) %>%
  tibble::rownames_to_column("kinship_consideration_group") %>%
  tidyr::pivot_longer(cols = starts_with("country"), values_to = "country") %>% 
  dplyr::select(-name) %>%
  dplyr::left_join(., dplyr::select(metadata, year, country, amplicon_barcode, genomics_sample_id))
SaveTabDelim(potential_links_specimens, paste0(output_dir, "/", batch_name, "_CountryBarcodeOverlapSpecimens.tsv"))


################################################################################
# All country pairwise similarity distributions
################################################################################
country_diff <- parallel::mclapply(setNames(barcode_countries, barcode_countries), function(i){
  diff <- SubsampleMerge(dplyr::filter(metadata, country == !!i & analysis_inclusion == "Included"))
}, mc.cores = 8) 
saveRDS(country_diff, file = paste0(output_dir, "/", batch_name, "_CountryPairwiseLists.rds"))

relatedness_summary <- lapply(country_diff, PairwiseWeightedSummary)
relatedness_df <- dplyr::bind_rows(lapply(relatedness_summary,"[[",1), .id="country_pair")  
SavePlots(similarity_distributions(relatedness_df), output_dir, "relatednessDensityMerge.png")


################################################################################
# Lab strain samples
################################################################################
lab_samples <- as.data.table(dplyr::filter(metadata, grepl("LAB", sample))) 
diff_lab <- SubsampleMerge(lab_samples, split_variable = NULL) 
SavePlots(related_and_missing_plot(diff_lab, "sample"), 
          output_dir, "lab_ClusterMatrices.png", width=6, height=5.5)



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
