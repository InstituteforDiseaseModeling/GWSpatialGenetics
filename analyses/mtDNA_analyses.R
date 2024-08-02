#!/usr/bin/env Rscript
################################################################################
# mtDNA Processing and Analysis Pipeline - A compilation of code
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: April 2024
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
  make_option(c("-p", "--site_missing"), type="numeric", default=0.75, 
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
packages_to_install <- c('vcfR', 'data.table', 'dplyr', 'tidyr', 'tibble', 
                         'stringr', 'geodist', 'ggplot2', 'ggpubr', 
                         'ggnewscale', 'ggridges')
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
  if(length(script.dir) == 0){
    print("Cannot determine script directory, setting script directory to working directory.")
    script.dir = getwd()
  }  
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

script_dir <- getScriptPath()
source(paste(script_dir, "module_checkMetadata.R", sep="/"))
source(paste(script_dir, "module_vcfdiploid2haploid.R", sep="/"))
source(paste(script_dir, "mtDNA_common_functions.R", sep="/"))


# plotitng options
theme_set(theme_bw())

# other functions
SaveTabDelim <- function(r_obj, path){
  write.table(r_obj, file = path, quote = F, row.names = F, sep="\t")
}


################################################################################
# input files and options
################################################################################
# set manual options for manual updates and runs of this code
diploid_vcf <- "/mnt/data/guinea_worm/processing/vcf_files/batch_July082024_jointGenotypeFiltered.vcf.gz"
metadata_file <- "/mnt/data/guinea_worm/metadata/GenomicSamplesMetadata_Database_v3.1_240611.txt"
project_dir <- "/mnt/data/guinea_worm/analyses"
# samp_missing <- 0.5
# site_missing <- 0.5
# min_gt_depth <- 5
# het_proportion  <- 1

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
vcf_clust <- gt2barcode(hap_vcf$vcf) %>%
  dplyr::mutate(vassar_worm = gsub("^[^_]*_|.batch.*", "", sample),
                filt_missing_prop = stringr::str_count(sequence, "\\.")/nchar(sequence),
                excluded_for_analysis = ifelse(filt_missing_prop > opt$filt_prop, "Yes", "No")) 

# get names of sampled in the vcf
vcf_samples <- colnames(hap_vcf$vcf@gt)[-1]
vcf_numbers <- cbind.data.frame(
  sample = vcf_samples,
  vassar_worm = gsub(".*_", "", vcf_samples))

# read in and check metadata
metadata <- data.table::fread(opt$metadata, sep="\t")
metadata[metadata == "." & !is.na(metadata)] <- NA
metadata <- run_metadata_checks(metadata, vcf_samples = vcf_numbers, excluded_samples = hap_vcf$missing_df)
metadata <- dplyr::left_join(metadata, vcf_clust) %>%
  dplyr::mutate(epi_foci = ifelse(epi_foci == "", NA, gsub(" focal area", "", epi_foci)),
                successfully_sequenced = ifelse(!is.na(filt_missing_prop), "Yes", "No"),
                host = factor(host, levels = names(host_colors)))


# save files for other analyses - i.e. Nextstrain, for collaborators
vcfR::write.vcf(hap_vcf$vcf,  paste0(output_dir, "/", batch_name, "_jointHaploidFilter.vcf.gz"))
SaveTabDelim(hap_vcf$summary, paste0(output_dir, "/", batch_name, "_filterSummary.tsv"))
SaveTabDelim(metadata, paste0(output_dir, "/", "GenomicSamplesMetadata_Database_", batch_name, ".tsv"))

# clear for memory
rm(vcf)

################################################################################
# One off requests
################################################################################
# flubendizole trial specimens - confirm they overlap and are D. medinensis 
# flu_samples <- read.delim("/mnt/data/guinea_worm/metadata/FLBZ_geneticsamples_2022.txt", header=F)
# flu_meta <- dplyr::filter(metadata, sample %in% flu_samples$V1) %>%
#   dplyr::select(sample, filt_missing_prop, excluded_for_analysis, successfully_sequenced)
# SaveTabDelim(flu_meta, paste0(output_dir, "/", "FlubendazoleSamples_", batch_name, ".tsv"))


################################################################################
# summaries
################################################################################
# missing positions to barcode frequency relationship 
barcode_missing <- dplyr::select(metadata, amplicon_barcode, sequence, frequency, filt_missing_prop) %>% unique()
missing_p <- barcode_missing %>%  ggplot(aes(x=filt_missing_prop, y=frequency)) +
  geom_jitter(alpha=0.5) +
  scale_y_log10() +
  geom_vline(xintercept = opt$filt_prop, linetype="dashed", color = "red") +
  labs(title = paste0(gsub("batch_", "Specimen library through ", batch_name), " sequencing batch"),       y = paste0("Frequency of barcode in sequenced library\n(", sum(barcode_missing$frequency, na.rm=T), " specimens)")) +
  annotate(geom = "text", color = "red",
           label = paste0("Sample missingness\nexclusion threshold\n(", sum(barcode_missing$filt_missing_prop > opt$filt_prop, na.rm=T), " specimens excluded)"),
           x = opt$filt_prop, y = 100, hjust = -0.1)
SavePlots(missing_p, output_dir, paste0(batch_name, "MissingFrequency.png"), 
          height=4, width=6)


# overall barcode distribution
metadata$amplicon <- factor(metadata$amplicon, levels=c(seq(1, 100), rev(c("Observed once", "Observed in < 5 samples", "Observed in < 10 samples", "Observed in < 20 samples"))))
specimen_counts <- table(metadata$country)
high_specimens_countries <- names(specimen_counts)[specimen_counts > 10]
barcode_counts <- metadata %>%
  dplyr::filter(!is.na(amplicon) & excluded_for_analysis == "No" & country %in% high_specimens_countries) %>%
  ggplot(aes(x=amplicon, y=..count.., fill=country)) +
  scale_fill_manual(values=country_colors) +
  geom_bar(color="black", size=0.25) +
  labs(x="Barcode group", y="Number of specimens", fill="Country") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
SavePlots(barcode_counts, output_dir, paste0(batch_name, "BarcodeHighCaseCountries.png"), 
          height=4, width=10)

# by country barcodes
barcode_countries <- names(table(metadata$country))[table(metadata$country) > 5]
lapply(barcode_countries, function(i){
  p_df <- metadata %>%
    dplyr::filter(country == !!i & excluded_for_analysis == "No") 
  
  # all 
  p <- p_df %>%
    ggplot(aes(x=year, y=..count.., fill=amplicon)) +
    geom_bar(color="black", size=0.25) +
    scale_fill_manual(values=nextstrain_colors) +
    scale_x_continuous(breaks = c(seq(min(p_df$year), max(p_df$year)))) +
    labs(title = paste("Barcodes in", i),
         subtitle = "(Groups determined by frequency in an all country library)",
         x="Year", y="Specimens", fill="Barcode") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  SavePlots(p, output_dir, paste0("barcode", gsub(" ", "", i), ".png"), 
            height=5, width=6)
  
  # proportion
  count_p <- p_df %>%
    ggplot(aes(x=year, y=..count.., fill=host)) +
    geom_bar(color="black", size=0.25) +
    scale_fill_manual(values=host_colors) +
    scale_x_continuous(breaks = c(seq(min(p_df$year), max(p_df$year)))) +
    labs(title = paste("Specimens from", i),
         x="Year", y="Specimens", fill="Host")
  prop_p <-p_df %>%
    ggplot(aes(x=year, y=..count.., fill=amplicon)) +
    geom_bar(position = "fill", color="black", size=0.25) +
    scale_fill_manual(values=nextstrain_colors) +
    scale_x_continuous(breaks = c(seq(min(p_df$year), max(p_df$year)))) +
    labs(title = paste("Barcodes in", i),
         subtitle = "(Groups determined by frequency in an all country library)",
         x="Year", y="Proportion of barcodes", fill="Barcode") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) 
  prop_comb <-  ggpubr::ggarrange(
    count_p + theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank()),
    prop_p, nrow=2, heights = c(.5,1), align = "v")
  SavePlots(prop_comb, output_dir, paste0("barcodeProportions", gsub(" ", "", i), ".png"), 
            height=8, width=6)
})

# barcodes by country overlap - groups mitochondrial lineages across long ranges 
count_barcodes <- dplyr::group_by(metadata, country, amplicon_barcode) %>%
  dplyr::summarise(country_count = n())
count_barcodes_total <- table(count_barcodes$amplicon_barcode)
shared_barcodes <- names(count_barcodes_total)[count_barcodes_total > 1]

potential_country_links <- dplyr::bind_rows(lapply(shared_barcodes, function(i){
  country_years <- dplyr::filter(metadata, amplicon_barcode == !!i) %>% 
    dplyr::group_by(country, year) %>% dplyr::summarise(country_yearly_specimens = n())
  # get years of each country in a list for quick comparison
  country_year_list <- with(dplyr::select(country_years, -country_yearly_specimens) %>% unique(), split(year, country))
  
  # loop to account for several countries and years per amplicon for kinship consideration
  cumulative_df <- data.frame()
  for(j in 1:length(country_year_list)){
    if(j < length(country_year_list)){
      country1 = names(country_year_list[j])
      country2 = names(country_year_list[j+1])
      overlapping_year = Reduce(intersect, list(a = country_year_list[[j]], 
                                                b = country_year_list[[j+1]]))
      if(length(overlapping_year > 1)){
        tmp <- cbind.data.frame(amplicon_barcode = as.integer(i), country1 = country1, country2 = country2, overlapping_year = overlapping_year) %>%
          dplyr::inner_join(., country_years %>% dplyr::rename("country1" = "country", "country1_yearly_specimens" = "country_yearly_specimens", "overlapping_year" = "year")) %>%
          dplyr::inner_join(., country_years %>% dplyr::rename("country2" = "country", "country2_yearly_specimens" = "country_yearly_specimens", "overlapping_year" = "year"))
        cumulative_df <- bind_rows(cumulative_df, tmp)
      } else{
        print(paste("No overlapping years for", country1, "and", country2, "for amplicon barcode", i))
      }
    }
  }
  return(cumulative_df)
})) %>% unique()

# save the overlapping barcodes between countries for the same year 
SaveTabDelim(potential_country_links, paste0(output_dir, "/", batch_name, "_CountryBarcodeOverlapCounts.tsv"))
# identify the specimens with shared mitochondrial DNA across countries for sibling prediction
potential_links_specimens <- dplyr::select(potential_country_links, amplicon_barcode, country1, country2, overlapping_year) %>%
  tibble::rownames_to_column("kinship_consideration_group") %>%
  tidyr::pivot_longer(cols = starts_with("country"), values_to = "country") %>% dplyr::select(-name) %>%
  dplyr::left_join(., dplyr::select(metadata, year, country, amplicon_barcode, genomics_sample_id))
SaveTabDelim(potential_links_specimens, paste0(output_dir, "/", batch_name, "_CountryBarcodeOverlapSpecimens.tsv"))


#######################################################################
# analysis data files - similarity
#######################################################################
# pairwise distance calculations
pairwise_dist <- reshape2::melt(HaversineDist(metadata)) %>% dplyr::rename(meters=value)

# genetic similarity calculations
tmp_gt <- data.frame(extract.gt(hap_vcf$vcf, element = "GT", return.alleles = TRUE))
gt <- data.frame(t(tmp_gt))
names(gt) <- row.names(tmp_gt)
row.names(gt) <- names(tmp_gt)
genetic_diff <- PairwiseDifference(gt)
genetic_diff <- dplyr::rename(genetic_diff, sample.x='Var1', sample.y='Var2')
SaveTabDelim(genetic_diff, paste0(output_dir, "/", batch_name, "_relatedness.txt"))

diff_df <- unique(dplyr::left_join(genetic_diff, pairwise_dist)) %>%
  mutate_cond(Var1==Var2, All=NA, rmNA=NA, missing=NA, meters=NA)
SaveTabDelim(diff_df, paste0(output_dir, "/", batch_name, "_relatedness.txt"))

diff_df <- fread(paste0(output_dir, "/", batch_name, "_relatedness.txt"))
metadata_slim <- dplyr::select(metadata, sample, host, country, year, emergence_date, wormnum_dpdxid) %>%
  dplyr::mutate(emergence_date = as.Date(emergence_date, "%m/%d/%y"))
diff_all <- diff_df %>%
  dplyr::rename(sample.x=Var1, sample.y=Var2) %>%
  dplyr::inner_join(., metadata_slim %>% dplyr::rename(sample.x=sample), by="sample.x") %>%
  dplyr::inner_join(., metadata_slim %>% dplyr::rename(sample.y=sample), by="sample.y") %>%
  dplyr::mutate(
    country_pair = ifelse(country.x == country.y, country.x, "Different countries"),
    country_pair = factor(country_pair, levels=names(country_colors)),
    host_pair = ifelse(host.x == host.y, host.x, "Different hosts"),
    year_pair = ifelse(year.x == year.y, year.x, "Across years"),
    diff_days = abs(as.numeric(difftime(emergence_date.x, emergence_date.y, unit="days"))),
    time_pair = ifelse(420 >= diff_days & diff_days >= 300, "10-14 Months", "Not in 10-14 months"))
rm(diff_df)

################################################################################
# analyses
################################################################################

################################################################################
# All country distributions
density_p <- similarity_distributions(diff_all %>% dplyr::filter(missing < opt$filt_prop & country_pair %in% c(barcode_countries, "Different countries")))
SavePlots(density_p, output_dir, "relatednessDensityMerge.png")
rm(density_p)


################################################################################
# Angola
ang <- dplyr::filter(metadata, country == "Angola")

ClusterPlots(dplyr::filter(ang, year > 2017), output_name="Angola2018+_PDB", sample_name = "wormnum_dpdxid")
SavePlots(dist_by_sim_plot(ang), output_dir, "simDist_ANG.png", height=5, width=10)

# Try new plot
devtools::install_github("jokergoo/circlize")
require("circlize")
diff_ang <- dplyr::filter(genetic_diff, str_detect(sample.x, "^ANG") & str_detect(sample.y, "^ANG")) %>% 
  MergeMeta(.)

chordDiagram(diff_ang %>%
  dplyr::rename(from=year.x, to=year.y) %>% 
  dplyr::select(from, to, rmNA), 
  link.visible = diff_ang$rmNA == 0)


# check which barcodes may be observed in other populations
angola_barcodes <- unique(ang$amplicon_barcode)
shared_barcodes <- dplyr::filter(metadata, metadata$amplicon_barcode %in% angola_barcodes)
unique(shared_barcodes$country)

# identify potential links for barcodes with missing positions
ang_missing <- dplyr::filter(ang, amplicon == "Observed once") %>% .[["sample"]]
ang_potential_links <- dplyr::filter(genetic_diff, (sample.x %in% ang_missing | sample.y %in% ang_missing) & rmNA == 0 & missing < opt$filt_prop) %>%
  MergeMeta(.)
if(nrow(dplyr::filter(ang_potential_links, country_pair != "Angola")) > 0){
  ClusterPlots(ang_potential_links, output_name="AngolaPotentialLinks", sample_name = "wormnum_dpdxid")
} else{
  print("No potential links for barcodes observed in current library.")
}


################################################################################
# Cameroon

# All samples
ClusterPlots(dplyr::filter(metadata, country == "Cameroon"), output_name="CameroonAll_Filter0.9+", filter_min=0.1)
ClusterPlots(dplyr::filter(metadata, country == "Cameroon" | adminb %in% c("Bongor", "Fianga")), filter_min=opt$filt_prop, output_name="Cameroon&BongorFianga_Filter0.9+")


# Cameroon and the northwestern Chadian border
border_sub <- dplyr::filter(metadata, country == "Cameroon" | adminb %in% c("Bongor", "Fianga"))

border_diff <- dplyr::filter(diff_all, sample.x %in% border_sub$sample & sample.y %in% border_sub$sample) %>%
  rowwise() %>%
  mutate(host_pair = ifelse(host_pair == " Different hosts", paste0(sort(c(host.x, host.y)), collapse = '-'), host_pair))

# counts 
count_p <- border_sub %>% ggplot(aes(x=as.character(year), y=..count.., fill=country))+
  geom_bar(alpha=0.5, color="black") +
  scale_fill_manual(values = country_colors) +
  labs(x="Year", y="Specimens", fill="Country")
SavePlots(count_p, output_dir, name="specimenCounts_Cameroon&BongorFianga.png", width=4, height=2)

# by range  
pairwise_counts_border <- border_diff %>% 
  dplyr::filter(missing < 0.1) %>%
  ggplot(aes(x=meters/1000,y=1-rmNA)) +
  #geom_hex() +
  geom_point(alpha=0.5) +
  labs(x="Distance between pairs of worm samples (kilometers)",
       y="Pairwise mtDNA genetic similarity\n(Excluding missing positions in barcode)",
       fill="Number of pairwise\ncomparisons") +
  facet_grid(~country_pair, scales = "free",
             labeller = labeller(country_pair = labels))
SavePlots(pairwise_counts_border, output_dir, "pairwiseHexByDist_Cameroon&BongorFianga.png", 
          height=3, width=10)

labels <- setNames(c("Within Cameroon", "Within Chad\n(Bongor and Fianga district)", 
                     "Cameroon and\nBongor or Fianga pair"), sort(unique(border_diff$country_pair))) 
density_border <- border_diff %>% 
  dplyr::filter(missing < 0.1) %>%
  ggplot(aes(x=1-rmNA, y=..scaled.., fill=country_pair)) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=country_colors,
                    labels = labels,
                    name="Country") +
  labs(x="Pairwise mtDNA genetic similarity\n(Excluding missing positions in barcode)", 
       y="Scaled probability density\nof pairwise comparisons")
SavePlots(density_border, output_dir, "relatednessDistribution_Cameroon&BongorFianga.png", 
          height=3, width=6)


################################################################################
# Ethiopia
eth <- dplyr::filter(metadata, country == "Ethiopia")
diff_eth <- dplyr::filter(diff_all, country_pair == "Ethiopia")

ClusterPlots(dplyr::filter(eth, year > 2012), output_name="Ethiopia2013+")
SavePlots(dist_by_sim_plot(eth), output_dir, "simDist_ETH.png", height=5, width=10)

eth_barcodes <- na.omit(unique(eth$amplicon_barcode))
eth_shared <- dplyr::filter(metadata, metadata$amplicon_barcode %in% eth_barcodes)

eth_shared_p <- eth_shared %>% ggplot(aes(x=year, y=..count.., fill=country))+
  geom_bar() +
  facet_wrap(~amplicon)
eth_shared_p

# Border regions
border_areas <- c("Abobo", "Gog")
ssu_samples <- c("SSUHUM2013_00495", "SSUHUM2013_00496", "SSUHUM2013_00507", "SSUHUM2021_04534")
ssu_border <- dplyr::filter(metadata, #zone %in% border_areas | 
                            genomics_sample_id %in% ssu_samples)

ssu_overlap <- dplyr::filter(metadata, amplicon %in% na.omit(unique(ssu_border$amplicon)))

# By host
diff_eth %>% ggplot(aes(x=1-rmNA)) + 
stat_ecdf(aes(color = host_pair)) +
  scale_color_manual(values=host_colors) +
  labs(x="Genetic similarity", y="Cumulative density",
       color="Country of pairs of worms")

# Specific baboon questions by Lexi 05/30/2024
sample_name <- "wormnum_dpdxid"
filter_min=0

# What is the relatedness is term of mitochondrial DNA for the Gog baboon worms, which is all the baboon from 2016-2020. 
# I have mapped the locations and presumed troops of each baboon. I know there were no sibling or parent worms between individual baboons,  but are there any “family” clusters that I might also compare to known movement and overlap between troops in this area? 
# How many mitochondrial types total have been found in baboons in Gog, and do any overlap with the baboons in Abobo (2022 baboons)?
eth_bab <- dplyr::filter(eth, host == "Baboon (Papio anubis)") 
df_cluster <- dplyr::filter(eth_bab, year > 2015 & year < 2023) 
gog_amp <- dplyr::filter(df_cluster, adminb == "Gog") %>% .[['amplicon_barcode']] %>% unique()
abobo_amp <- dplyr::filter(df_cluster, adminb == "Abobo") %>% .[['amplicon_barcode']] %>% unique()
intersect(gog_amp, abobo_amp)

tmp_pair <- dplyr::filter(diff_eth, sample.x %in% df_cluster$sample & sample.y %in% df_cluster$sample) %>%
  mutate(outline = ifelse(rmNA == 0, T, NA))
related_joint <- related_and_missing_plot(tmp_pair, eval(sample_name), filter_min=filter_min)
SavePlots(related_joint, output_dir, "ETHGogandAboboRelatedness.png", 
          height=15, width=15)

bab_amps <- gog_bab %>%
  dplyr::mutate(amplicon2 = ifelse(frequency == 1, "Observed once", amplicon_barcode)) %>%
  ggplot(aes(x=year, y=..count.., fill=amplicon2)) +
  geom_bar(color="black", size=0.25) +
  #scale_fill_manual(values=nextstrain_colors) +
  scale_x_continuous(breaks = c(seq(min(eth_bab$year), max(eth_bab$year)))) +
  labs(title = paste("Barcodes in Ethiopian Baboons (Gog and Abobo)"),
       subtitle = "(Groups determined by frequency in an all country library)",
       x="Year", y="Specimens", fill="Barcode") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
SavePlots(bab_amps, output_dir, "ETHGogandAboboAmplicons.png", 
          height=6, width=6)


# What is the relatedness is term of mitochondrial DNA for Baby (dog, 2020), the 2 Gutok baboons (2022), the leopard (2019), and the serval (2023)? I know most have multiple non-sibling worms. 
# All these worms were found in a small area (within ~20 km), and there is likely only one or two troops between Baby and the serval from known baboon movement.  
baby <- dplyr::filter(eth, animal_name == "Baby") %>% .[['genomics_sample_id']]
multiple_worm_ids <- dplyr::filter(eth, grepl("ETHBAB2022|ETHPPD2019|ETHWCT2023", genomics_sample_id)) %>% .[['genomics_sample_id']]
filter_ids <- c(baby, multiple_worm_ids)
df_cluster <- dplyr::filter(eth, genomics_sample_id %in% filter_ids)
tmp_pair <- dplyr::filter(diff_eth, sample.x %in% df_cluster$sample & sample.y %in% df_cluster$sample) %>%
  mutate(outline = ifelse(rmNA == 0, T, NA))
related_joint <- related_and_missing_plot(tmp_pair, eval(sample_name), filter_min=filter_min)
SavePlots(related_joint, output_dir, "ETHSmallAreaRelatedness.png", 
          height=8, width=8)



################################################################################
# South Sudan
ssu <- dplyr::filter(metadata, country == "South Sudan")
diff_ssu <- dplyr::filter(diff_all, country_pair == "South Sudan")
# add in epi focus
diff_ssu <-  diff_ssu %>% 
  dplyr::inner_join(., dplyr::select(metadata, sample, epi_foci) %>% dplyr::rename(sample.x=sample), by="sample.x") %>%
  dplyr::inner_join(., dplyr::select(metadata, sample, epi_foci) %>% dplyr::rename(sample.y=sample), by="sample.y") %>%
  rowwise() %>%
  dplyr::mutate(focal_area = ifelse(is.na(epi_foci.x) | is.na(epi_foci.y), "Unknown foci for one\nor both worms in a pair",
                                    ifelse(epi_foci.x == epi_foci.y, epi_foci.x, paste0(sort(c(epi_foci.x, epi_foci.y)), collapse = '-'))),
                focal_area_2 = ifelse(!grepl("-", focal_area), "Same focal area", "Between focal areas"),
                focal_area_3 = ifelse(focal_area_2 == "Between focal areas", focal_area_2, focal_area),
                year_pair_2 = ifelse(year_pair == "Across years" & abs(year.x-year.y) == 1, paste0(sort(c(as.character(year.x), as.character(year.y))), collapse = '-'), year_pair),
                genetic_pair = ifelse(rmNA == 0, "Potentially linked specimens", "Likely unlinked specimens"),
                genetic_pair = ifelse(missing < 0.1, genetic_pair, "Undetermined due to worm\npair sequencing quality"),
                genetic_pair_2 = ifelse(year_pair != "Across years", paste0(genetic_pair, "- Same year"), paste(genetic_pair, "- Different year")))


# pairwise similarity plots
ClusterPlots(dplyr::filter(ssu, year > 2017), output_name="SouthSudan2018+_PDB", sample_name = "wormnum_dpdxid")
ClusterPlots(dplyr::filter(ssu, year == 2022), output_name="SouthSudan2022_PDB", sample_name = "wormnum_dpdxid")

# Focal comparisons
ssu_foci <- unique(ssu$epi_foci)
lapply(ssu_foci[-length(ssu_foci)], function(x){
  ClusterPlots(dplyr::filter(ssu, epi_foci == x), output_name=paste("SouthSudan_foci", gsub(" ", "", x)))#, sample_name = "wormnum_dpdxid")
})

ssu_foci_counts <- ssu %>% ggplot(aes(x=year, y=..count.., fill=epi_foci)) +
  geom_bar() +
  stat_count(geom = "text", colour = "black", size = 2.5,
             aes(label = ..count..),position=position_stack(vjust=0.5)) +
  scale_x_continuous(breaks = ~round(unique(pretty(.)))) +
  labs(x="Year", y="Number of specimens", fill="Region") 

diff_by_year <- dplyr::filter(diff_ssu, sample.x != sample.y & focal_area %in% ssu_foci & year_pair_2 != "Across years") %>%
  ggplot(aes(x=year_pair_2, y=..count.., fill=genetic_pair)) +
  geom_bar() +
  facet_grid(focal_area~., scales="free_y") +
  labs(x="Year of infection between pairs of worms", y="Pairwise comparisons", fill="")
diff_by_year

ssu_order <- c("Northern Jonglei", "East of Nile", "Central", "West of Nile", 
               "Between focal areas", "Unknown foci for one\nor both worms in a pair") 
diff_compare_dis <- dplyr::filter(diff_ssu, sample.x != sample.y & 
                                year_pair_2 != "Across years" & !grepl(year_pair_2, "-")) %>%
  dplyr::mutate(focal_area_3 = factor(focal_area_3, levels=ssu_order)) %>%
  ggplot(aes(x=focal_area_3, y=..count.., fill=genetic_pair)) +
  geom_bar() +
  labs(x="Focal area", y="Pairwise comparisons", fill="")

diff_compare_cont <- dplyr::filter(diff_ssu, (sample.x != sample.y) & 
                                     (year_pair_2 != "Across years") & (!grepl(year_pair_2, "-") & missing < opt$filt_prop)) %>%
  dplyr::mutate(focal_area_3 = factor(focal_area_3, levels=ssu_order),
                year_pair_3 = ifelse(year_pair != "Across years", "Same year", "Different years"),
                year_pair_3 = factor(year_pair_3, levels = c("Same year", "Different years"))) %>%
  ggplot(aes(y=focal_area_3, x=1-rmNA, linetype=year_pair_3)) +
  ggridges::geom_density_ridges(quantile_lines = TRUE, quantiles = 2, scale = 0.9, alpha=0.1) +
  coord_flip() +
  xlim(0.85,1) +
  labs(y="Focal area", x="Pairwise genetic similarity\n(Excluding missing positions)", linetype="Emergence year between\npairs of worms")
ssu_foci_merge <- ggpubr::ggarrange(diff_compare_dis, diff_compare_cont, ncol=1, align = "v")
SavePlots(ssu_foci_merge, output_dir, "focalComparison_SSU.png", height=6, width=9)


# dog sample comparison
ssu_dog <- dplyr::filter(ssu, host == "Dog")
ClusterPlots(ssu_dog, output_name="SouthSudanDogs", sample_name = "wormnum_dpdxid")

for (i in unique(ssu_dog$sample)){
  tmp_dog <- dplyr::filter(diff_ssu, sample.x == i | sample.y==i) %>%
    dplyr::filter(rmNA ==0 )
  subset_samples <- unique(c(tmp_dog$sample.x, tmp_dog$sample.y))
  tmp_df <- dplyr::filter(ssu, sample %in% subset_samples)
  # check for duplicates - happens with samples that did not have exact matching
  dups <- tmp_df %>% 
    group_by(wormnum_dpdxid) %>% 
    mutate(num_dups = n())
  if(any(dups$num_dups > 1)){
    duplicates <- dplyr::filter(dups, num_dups >1) %>% .[["wormnum_dpdxid"]] %>% unique()
    dup_df <- data.frame()
    for(dup in duplicates){
      tmp_dup <- dplyr::filter(tmp_df, wormnum_dpdxid == dup) 
      tmp_dup <- tmp_dup[which.max(tmp_dup$amplicon_barcode),]
      dup_df <- rbind(dup_df, tmp_dup)
    }
    tmp_df <- bind_rows(dplyr::filter(tmp_df, !wormnum_dpdxid %in% duplicates), dup_df)
  }
  ClusterPlots(tmp_df, output_name=i, sample_name = "wormnum_dpdxid")
}

dist_p <- dist_by_sim_plot(ssu)
SavePlots(dist_p, output_dir, "simDist_SSU.png", height=5, width=10)


################################################################################
# Chad
# This will require updates as the genomics working group comes up with a more detailed analytical plan 
chad <- dplyr::filter(metadata, country == "Chad")
chad_admins_counts <- table(chad$admina)
chad_admins <- names(chad_admins_counts)[chad_admins_counts > 10]

SavePlots(dist_by_sim_plot(chad), output_dir, "simDist_CHD.png", height=5, width=10)

lapply(chad_admins, function(i){
  p_df <- metadata %>%
    dplyr::filter(admina == !!i & excluded_for_analysis == "No") %>%
  # dplyr::filter(adminb %in% c("Bongor", "Fianga") & excluded_for_analysis == "No" & !is.na(amplicon)) %>%
  # p_df <- dplyr::filter(chad, admina %in% chad_admins & !is.na(amplicon)) %>%
    dplyr::mutate(amplicon = ifelse(frequency > 11 & frequency < 20, "Observed in < 20 samples", as.character(amplicon)))
  p_df$amplicon <- factor(p_df$amplicon, levels=c(seq(1, 100), rev(c("Observed once", "Observed in < 5 samples", "Observed in < 10 samples", "Observed in < 20 samples"))))
  
  # all 
  if(nrow(p_df) > 1){
   p <- p_df %>%
    ggplot(aes(x=year, y=..count.., fill=amplicon)) +
    geom_bar(color="black", size=0.25) +
    scale_fill_manual(values=nextstrain_colors) +
    scale_x_continuous(limits = c(2012, 2024), breaks = seq(2012, 2024)) +
    labs(title = paste("Barcodes in", i),
         subtitle = "(Groups determined by frequency in an all country library)",
         x="Year", y="Specimens", fill="Barcode") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    facet_wrap(~admina, scales="free_y")
  SavePlots(p, output_dir, paste0("barcodeChadAdmin", gsub(" ", "", i), ".png"), height=5, width=6)
  # SavePlots(p, output_dir, paste0("barcodeChadAdminFacet.png"), height=5, width=10)
  }
})
