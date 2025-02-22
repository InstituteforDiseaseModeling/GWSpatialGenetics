################################################################################
# mtDNA Processing and Analysis Pipeline 
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: February 2025
# Requires inputs created from running module_mtDNAAnalyeses.R to generate barcodes (mitochondrial lineage groupings) and pairwise barcode calculations. 
################################################################################

################################################################################
# input files and options
################################################################################
# set command line and default options
require('optparse')
option_list = list(
  make_option(c("-m", "--metadata"), type="character", 
              default=NULL, 
              help="Metadata processed from barcode analyses.", metavar="character"),
  make_option(c("-p", "--country_pairwise"), type="character", 
              default=NULL, 
              help="R object that contains all relevant pairwise_comparisone", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=paste0(format(Sys.time(), "%Y%m%d"), "_run"), 
              help="Output directory. A subdirectory with the batch name pulled from the  If none provided, will output to a folder in the running directory with the date.", metavar="character"),
  make_option(c("-f", "--filt_prop"), type="numeric", default=0.1, 
              help="The proportion of missing positions, after filtering high quality variants, tolerated in a barcode for analyses. Default is 10% of the barcode length.", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


################################################################################
# set-up
################################################################################
# load libraries
packages_to_install <- c('data.table', 'dplyr', 'tidyr', 'tibble', 'stringr', 
                         'ggplot2', 'ggpubr','ggnewscale', 'ggridges')
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

#script_dir <- getScriptPath()
setwd('/home/jribado/git/GWSpatialGenetics/analyses')
script_dir <- try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
source(paste(script_dir, "mtDNA_common_functions.R", sep="/"))
source(paste(script_dir, "mtDNA_plotting.R", sep="/"))

# other functions
reformat_cluster <- function(cluster_df, sheet, metadata_df = metadata){
  tmp <- readxl::read_excel(cluster_df, sheet = sheet) %>%
    dplyr::mutate(vassar_worm = gsub(".*_", "", `Genomics Sample ID`)) %>%
    dplyr::left_join(., dplyr::select(metadata_df, sample, amplicon_barcode, amplicon, frequency, wormnum_dpdxid, vassar_worm, host, frequency))
}



ClusterPariwiseMatrices <- function(cluster_name, df = metadata){
  cluster_df <- dplyr::filter(df, cluster == !!cluster_name & analysis_inclusion == "Included") 
  if(nrow(cluster_df) > 2){
    cluster_diff <- SubsampleMerge(cluster_df, split_variable = NULL)
    ClusterPlots(df = cluster_df, pairwise_df = cluster_diff,
                 output_name=paste0("cluster_", gsub(" |cluster", "", cluster_name)))
  }
}


# plotting options
theme_set(theme_bw())


################################################################################
# Load data and set variables
################################################################################
# metadata_file = '/mnt/data/guinea_worm/analyses/batch_Dec052024/GenomicSamplesMetadata_Database_batch_Dec052024.rds'

# update parser if running in manual mode
if(is.null(opt$metadata)){
  opt$metadata <- metadata_file
  opt$output_dir <- "/mnt/data/guinea_worm/analyses"
}

batch_name <- gsub(".*Database_|*\\.rds", "", basename(opt$metadata))
output_dir <- paste(opt$output_dir, batch_name, sep="/")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

if(grepl("\\.txt$", opt$metadata)){
  metadata <- data.table::fread(opt$metadata, sep="\t")
} else{
  metadata <- readRDS(opt$metadata)
} 

country_diff <- readRDS(paste0(dirname(opt$metadata), "/", batch_name, "_CountryPairwiseLists.rds"))
barcode_pairwise <- data.table::fread(paste0(dirname(opt$metadata), "/", batch_name, "_relatednessBarcodes.txt"))
missing_count_variable = names(metadata)[grepl("missing_n_", names(metadata))]
n_variants <- readr::parse_number(missing_count_variable)

# update color for amplicons 
common_amplicons <- readr::parse_number(
  as.character(unique(metadata$amplicon)[!unique(metadata$amplicon) %in% names(reduced_base)])) 
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

################################################################################
# Barcode rarefaction 
# In progress 02/2025
# Goals: Provide an upper bound of missed lineages per year to use as baseline for expected lineages the following year
################################################################################
inext_res <- LineageExtrapolation("Chad")

################################################################################
# Epidemiological vignettes
# Porting from 2023 PRM analyses - may need updating 
# These vignettes have been informed in collaboration with TCC epidemiologists Maryann Delea and Obiora Eneanya. 
# * Chad
#   - Bogam cluster (specimens from 2019-2021) 
# * Ethiopia
#   - Duli Farm cluster (baboon and human specimens from 2019-2021) 
#   - 2018 PRC Agnuak cat cluster 
#   - 2020 PRC Aguak cat cluster 
#   - 2020 Ogul Ponds cluster
################################################################################
old_clusters <- dplyr::bind_rows(
  reformat_cluster("/mnt/data/guinea_worm/metadata/GW_clusters_genomics_30Jan2023.xlsx", 3)%>%
    dplyr::mutate(collective = as.character(collective)),
  reformat_cluster("/mnt/data/guinea_worm/metadata/GW_clusters_genomics_30Jan2023.xlsx", 1) %>%
  dplyr::mutate(cluster = gsub("Link to |\\ cluster", "", cluster))
) %>%
  dplyr::mutate(cluster = ifelse(cluster %in% c("88", "99"), NA, cluster))

metadata <- dplyr::left_join(metadata, dplyr::select(old_clusters, sample, cluster))

lapply(unique(metadata$cluster)[!is.na(unique(metadata$cluster))], ClusterPariwiseMatrices)


################################################################################
# Angola
################################################################################
ClusterPlots(dplyr::filter(dplyr::filter(metadata, country == "Angola"), year > 2022), 
             pairwise_df = country_diff['Angola'][[1]], 
             output_name="Angola2022+_PDB", sample_name = "wormnum_dpdxid")
SavePlots(dist_by_sim_plot(dplyr::filter(metadata, country == "Angola"), 
          country_diff["Angola"][[1]]), 
          output_dir, "ANG_simDist.png", height=5, width=10)

relatedness_summary <- PairwiseWeightedSummary(country_diff["Angola"][[1]], by_year = T)
SavePlots(PairwiseProportionHeatMap(relatedness_summary$weighted_diff),
          output_dir, "ANG_pairwiseHeatmap.png", height=5, width=5)

# check which barcodes may be observed in other populations
angola_barcodes <- dplyr::filter(metadata, country == "Angola") %>% dplyr::distinct(amplicon_barcode)
shared_barcodes <- dplyr::filter(metadata, metadata$amplicon_barcode %in% angola_barcodes)
if(length(unique(shared_barcodes$country)) > 1){
  print("Countries with an overlapping Angola barcode:", unique(shared_barcodes$country)[unique(shared_barcodes$country) != "Angola"])
}

ang_potential_links <- dplyr::filter(barcode_pairwise, 
  (amplicon_barcode.x %in% angola_barcodes | amplicon_barcode.y %in% angola_barcodes) & 
   amplicon_barcode.x != amplicon_barcode.y &  
   rmNA == 0 & missing < opt$filt_prop)
if(nrow(ang_potential_links) > 0){
  print("Some potential outside Angola links identified.")
  ang_shared_barcode_samples <- dplyr::filter(metadata, amplicon_barcode %in% angola_barcodes & country != "Angola")
} 


################################################################################
# Cameroon
################################################################################  
# All samples
relatedness_summary <- PairwiseWeightedSummary(country_diff["Cameroon"][[1]], by_year=T)
SavePlots(PairwiseProportionHeatMap(relatedness_summary$weighted_diff),
          output_dir, "CAM_pairwiseHeatmap.png", height=5, width=5)


ClusterPlots(dplyr::filter(metadata, country == "Cameroon"), 
             pairwise_df = country_diff['Cameroon'][[1]],
             output_name="CameroonAll_Filter0.9+", filter_min=0.1)


# Cameroon and the northwestern Chadian border
####################################
border_sub  <- dplyr::filter(metadata, country == "Cameroon" | adminb %in% c("Bongor", "Fianga"))
# border_sub <- dplyr::filter(metadata, country == "Cameroon" | epi_foci == "Chad-Cameroon border area")
diff_border <-SubsampleMerge(border_sub, split_variable = NULL)

cam_border_distribution <- VariableDistributionPlot(
  df = border_sub, dt = diff_border, 
  variable= "country", color_vector = country_colors) 
cam_border_distribution <- ggpubr::annotate_figure(cam_border_distribution, 
  top = "Pairwise relatedness in the Cameroon and\nChad border ( (Bongor and Fianga district)")
SavePlots(cam_border_distribution, output_dir, "CAM-CHD_borderPairwiseDistributions.png", height=4, width=6)


# by range  
labels <- setNames(c("Within Cameroon", "Within Chad\n(Bongor and Fianga district)", 
                     "Cameroon and\nBongor or Fianga pair"), sort(unique(diff_border$country_pair))) 
pairwise_counts_border <- diff_border %>% 
  dplyr::filter(missing < 0.1) %>%
  ggplot(aes(x=meters/1000,y=1-rmNA)) +
  #geom_hex() +
  geom_point(alpha=0.5) +
  labs(x="Distance between pairs of worm samples (kilometers)",
       y="Pairwise mtDNA genetic similarity\n(Excluding missing positions in barcode)",
       fill="Number of pairwise\ncomparisons") +
  facet_grid(~country_pair, scales = "free",
             labeller = labeller(country_pair = labels))
SavePlots(pairwise_counts_border, output_dir, "CAM-CHD_pairwiseHexByDist.png", 
          height=3, width=8)


################################################################################
# Ethiopia
################################################################################
SavePlots(dist_by_sim_plot(
  dplyr::filter(metadata, country == "Ethiopia"), 
  country_diff['Ethiopia'][[1]]), 
  output_dir, "ETH_simDist.png", height=5, width=10)

relatedness_summary <- PairwiseWeightedSummary(country_diff["Ethiopia"][[1]], by_year = T)
SavePlots(PairwiseProportionHeatMap(relatedness_summary$weighted_diff),
          output_dir, "ETH_pairwiseHeatmap.png", height=5, width=5)


# Border regions
####################################
border_areas <- c("Abobo", "Gog")
ssu_samples <- c("SSUHUM2013_00495", "SSUHUM2013_00496", "SSUHUM2013_00507", "SSUHUM2021_04534")
ssu_border <- dplyr::filter(metadata, adminb %in% border_areas | 
                            genomics_sample_id %in% ssu_samples)


# Host comparisons
####################################
eth_host_distribution <- VariableDistributionPlot(
  df = dplyr::filter(metadata, country == "Ethiopia"),
  dt = country_diff['Ethiopia'][[1]], 
  variable= "host", color_vector = host_colors) 
eth_host_distribution <- ggpubr::annotate_figure(eth_host_distribution, 
  top = "Pairwise relatedness of infected host in Ethiopia")
SavePlots(eth_host_distribution, plot_path = output_dir, name = "ETH_hostPairwiseDistributions.png", height = 4, width=8)

# host barcodes
BarcodeCountPlots(dplyr::filter(metadata, country == "Ethiopia" & host %in% c("Human", "Baboon (Papio anubis)") & !is.na(amplicon)),
                  ggplot_title = "ETH Host", facet_variable = "host")

eth_baboon <- dplyr::filter(metadata, country == "Ethiopia" & host %in% c("Human", "Baboon (Papio anubis)") & !is.na(amplicon)) %>%
  separate(wormnum_dpdxid, c('PD_ID', 'Worm'), sep="\\.", remove = F)
eth_bab_counts <- sort(table(eth_baboon$PD_ID), decreasing = T)
baboon_barcodes_p <- eth_baboon %>%
  dplyr::filter(PD_ID %in% names(eth_bab_counts[eth_bab_counts > 1])) %>%
  dplyr::mutate(PD_ID = factor(PD_ID, levels = names(eth_bab_counts)))  %>%
  ggplot(aes(x=PD_ID, y=after_stat(count), fill=amplicon, color=as.character(amplicon_barcode))) +
  geom_bar(color="black") +
  scale_fill_manual(values=nextstrain_colors) +
  labs(x=paste0("Unique baboons\n(Program PD, ", eth_bab_counts[eth_bab_counts > 1],
                " baboons with more than one specimen from\n", length(unique(eth_baboon$PD_ID)), " total baboons)"),
       y=paste0("Worm specimens\n(n=", nrow(eth_baboon), " total)"),
       fill="Barcode") +
  guides(x =  guide_axis(angle = 45))
SavePlots(baboon_barcodes_p , plot_path = output_dir, name = "ETH_BaboonBarcodes.png", height = 4, width=6.5)

ClusterPlots(df = dplyr::filter(metadata, country == "Ethiopia" & host == "Baboon (Papio anubis)"),
            pairwise_df = SubsampleMerge(dplyr::filter(metadata, country == "Ethiopia" & host == "Baboon (Papio anubis)"), split_variable = NULL),
            output_name =  "ETH_BaboonRelatednessMatrix")


################################################################################
# South Sudan
################################################################################
# Overall pairwise distance by relatedness relationship
####################################
SavePlots(dist_by_sim_plot(
  dplyr::filter(metadata, country == "South Sudan"), 
  country_diff['South Sudan'][[1]]), 
  output_dir, "SSU_simDist.png", height=5, width=10)

relatedness_summary <- PairwiseWeightedSummary(country_diff["South Sudan"][[1]], by_year = T)
SavePlots(PairwiseProportionHeatMap(relatedness_summary$weighted_diff),
          output_dir, "SSU_pairwiseHeatmap.png", height=5, width=5)


# Focal comparisons
####################################
BarcodeCountPlots(dplyr::filter(metadata, country == "South Sudan") %>%
                    dplyr::mutate(epi_foci = ifelse(is.na(epi_foci), "Not provided", epi_foci),
                                  epi_foci = factor(epi_foci, levels=c(names(foci_order), "Not provided"))),
                  facet_variable = "epi_foci")

ssu_focal_distribution <- VariableDistributionPlot(
    df = dplyr::filter(metadata, country == "South Sudan"),
    dt = country_diff['South Sudan'][[1]], 
    variable= "epi_foci") 
ssu_focal_distribution <- ggpubr::annotate_figure(ssu_focal_distribution, 
    top = "Pairwise relatedness in South Sudan by focal areas")
SavePlots(ssu_focal_distribution, output_dir, "SSU_focalPairwiseDistributions.png", height=3, width=7)


# Dog specimen comparisons
####################################
ssu_dog <- dplyr::filter(metadata, country == "South Sudan" & host == "Dog") %>% dplyr::pull(amplicon_barcode) %>% unique()
SaveTabDelim(dplyr::filter(metadata, amplicon_barcode %in% ssu_dog), paste(output_dir, "SSU_DOGLineagesShared.txt", sep="/"))
             

################################################################################
# Chad
# This will require updates as the genomics working group comes up with a more detailed analytical plan 
SavePlots(dist_by_sim_plot(dplyr::filter(metadata, country == "Chad"),
                                         country_diff['Chad'][[1]]), 
          output_dir, "CHD_simDist.png", height=5, width=10)

relatedness_summary <- PairwiseWeightedSummary(country_diff["Chad"][[1]], by_year = T)
SavePlots(PairwiseProportionHeatMap(relatedness_summary$weighted_diff),
          output_dir, "CHD_pairwiseHeatmap.png", height=6, width=10)

# Admin specific comparison
chad_admins_counts <- dplyr::filter(metadata, country == "Chad") %>% pull(admina) %>% table()
chad_admins <- names(chad_admins_counts)[chad_admins_counts > 10]
lapply(chad_admins, function(i){
  BarcodeCountPlots(dplyr::filter(metadata, admina == !!i & analysis_inclusion == "Included"),
                    ggplot_title = paste0("Chad -", i))
})  
  
chd_admin_distribution <- VariableDistributionPlot(
  df = dplyr::filter(metadata, country == "Chad"),
  dt = country_diff['Chad'][[1]], 
  variable= "admina") 
chd_admin_distribution <- ggpubr::annotate_figure(chd_admin_distribution, top = "Pairwise relatedness of regions in Chad")
SavePlots(chd_admin_distribution, plot_path = output_dir, name = "CHD_adminaPairwiseDistributions.png", height = 4, width=8)
 
