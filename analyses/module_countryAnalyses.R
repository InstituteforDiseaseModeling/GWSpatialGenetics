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
                         'ggplot2', 'ggpubr','ggnewscale', 'ggdist', 'ggridges', 'magrittr')
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
# getScriptPath <- function(){
#   cmd.args <- commandArgs()
#   m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
#   script.dir <- dirname(regmatches(cmd.args, m))
#   if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
#   if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
#   return(script.dir)
# }
# 
# script_dir <- getScriptPath()
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


FirstBarcodeObservation <- function(input_key){
  barcodes_by_year  <- metadata[admin_key == input_key, ..columns] %>% group_by_at(columns) %>% summarize(n_specimens = n())
  barcode_vector <- split(barcodes_by_year$year, barcodes_by_year$amplicon_barcode)
  barcodes_by_year <- rowwise(barcodes_by_year) %>%
    dplyr::mutate(barcode_emergence = case_when(
      (as.numeric(year)-1) %in% barcode_vector[match(as.character(amplicon_barcode), names(barcode_vector))][[1]] ~ "Seen in previous year",
      any(barcode_vector[match(as.character(amplicon_barcode), names(barcode_vector))][[1]] < as.numeric(year) - 1 ) ~ "Seen in > 1 previous years",
      amplicon_barcode %in% singletons ~ "Observed once, cannot be linked", 
      TRUE ~ "First barcode emergence"))
  return(barcodes_by_year)   
}


NearbyBarcodeRanges <- function(input_key, data = admin_barcode_emergence){
  data <- dplyr::filter(data, admin_key == input_key & barcode_emergence == "First barcode emergence")
  first_year_vector <- setNames(data$year, data$amplicon_barcode)
  admin_gps_ids <- unique(metadata[admin_key == input_key, gps_id])
  
  dist_summary <- lapply(names(first_year_vector), function(i){
    barcode_sub <- metadata[country == "Chad" & amplicon_barcode == as.numeric(i) & year == first_year_vector[i][[1]] - 1]
    if(nrow(barcode_sub) > 1){
      amplicon_pairs <- expand.grid(unique(barcode_sub$gps_id), admin_gps_ids) %>%
        dplyr::mutate(gps_id_pair = paste(pmin(Var1, Var2), pmax(Var1, Var2), sep="_"),
                      Var1 = as.character(Var1)) %>%
        dplyr::inner_join(., as.data.frame(table(barcode_sub$gps_id), stringsAsFactors=F)) %>%
        dplyr::left_join(., chad_gps_pairs)
      
      amplicon_summary <- amplicon_pairs %>%
        summarize(barcode_specimens = sum(Freq),
                  minimum  = min(meters, na.rm = T),
                  maximum  = max(meters, na.rm = T),
                  weighted_mean = weighted.mean(meters, Freq, na.rm = T),
                  weighted_median = weighted.median(meters, Freq))
    
      amplicon_summary$admin_key = input_key
      amplicon_summary$amplicon_barcode = i
      amplicon_summary$previous_year = first_year_vector[i] - 1
      return(amplicon_summary)
    }
  })
  if(length(dist_summary) > 1){
    admin_summary <- bind_rows(dist_summary)
    return(admin_summary)
  }
}


GroupByCounts <- function(metadata, amplicon_group = "amplicon_barcode", subset_country = "Chad", by_village=FALSE){
  subset_columns <- c("country", "year", amplicon_group, missing_count_variable, "analysis_inclusion")
  if(isTRUE(by_village)){
    admin_columns  <- c("admina", "adminb", "adminc", "admind")
    subset_columns <- c(subset_columns, admin_columns)
  }
  
  if(amplicon_group != "amplicon_barcode_conservative"){
    group_by_columns <- subset_columns[subset_columns != missing_count_variable]
    missing_df <- dplyr::select(metadata, amplicon_barcode_conservative, eval(missing_count_variable)) %>% unique() %>% 
      dplyr::group_by(amplicon_barcode_conservative) %>% 
      dplyr::filter(eval(missing_count_variable) == min(eval(missing_count_variable))) %>%
      dplyr::rename("amplicon_barcode"="amplicon_barcode_conservative")
  } else{
    group_by_columns <- subset_columns
  }
  
  tmp_counts <- dplyr::select_(metadata, .dots = subset_columns) %>% 
    dplyr::filter(country == !!subset_country) %>%
    dplyr::group_by_at(group_by_columns) %>%
    summarise(n_specimens=n(), .groups = 'drop')
  
  if(amplicon_group != "amplicon_barcode_conservative"){
    tmp_counts <- dplyr::left_join(tmp_counts, missing_df)
  }
  return(tmp_counts)
}

# plotting options
theme_set(theme_bw())


################################################################################
# Load data and set variables
################################################################################
# metadata_file = '/mnt/data/guinea_worm/analyses/batch_Feb242025/GenomicSamplesMetadata_Database_batch_Feb242025.rds'

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

metadata <- dplyr::left_join(unique(metadata), dplyr::select(old_clusters, sample, cluster))

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
diff_border = SubsampleMerge(border_sub, split_variable = NULL)

BarcodeCountPlots(border_sub, facet_variable = "country", 
                  ggplot_title = "CAM-CHD", input_width=8, input_height=4.5)

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

# Work in progress - DAPC analyses
# vcf <- vcfR::read.vcfR(paste0(output_dir, "/", batch_name, "_jointHaploidFilterWithSingletons.vcf.gz"))
# border_vcf <- vcf[sample=border_sub$sample]
# border_sub <- border_sub[match(colnames(border_vcf@gt)[-1], border_sub$sample),]
# my_genind <- vcfR::vcfR2genind(border_vcf)
# my_genind@pop <- as.factor(border_sub$country)
# 
# grp <- find.clusters(my_genind, max.n.clust=40)
# table(pop(my_genind), grp$grp)
# dapc1 <- dapc(my_genind, pop=pop(my_genind), n.pca=20, n.da=6)
# ggsave(filename=paste(output_dir, "CAM-CHD_DAPCByPop.png", sep="/"),
#   plot = call(scatter(dapc1, scree.da=FALSE,  legend=TRUE, col=c("red", "blue"), solid=.4)), 
#   device = "png")
# 
# 
# compoplot(dapc1, posi="bottomright",
#           txt.leg=paste("Cluster", 1:2), lab="",
#           ncol=1, xlab="individuals", col=funky(2),
#           )

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

eth_baboon <- dplyr::filter(metadata, country == "Ethiopia" & grepl("BAB", sample)) %>%
  separate(wormnum_dpdxid, c('PD_ID', 'Worm'), sep="\\.", remove = F)
eth_bab_counts <- sort(table(eth_baboon$PD_ID), decreasing = T)
baboon_barcodes_p <- eth_baboon %>%
  dplyr::filter(PD_ID %in% names(eth_bab_counts[eth_bab_counts > 1])) %>%
  dplyr::mutate(PD_ID = factor(PD_ID, levels = names(eth_bab_counts)))  %>%
  ggplot(aes(x=PD_ID, y=after_stat(count), fill=amplicon, color=amplicon_barcode, group=amplicon_barcode)) +
  geom_bar(color = "black", linewidth = 0.5) +  # Adjust outline thickness as needed
  scale_fill_manual(values = nextstrain_colors) +
  labs(x=paste0("Unique baboons\n(Program PD, ", length(eth_bab_counts[eth_bab_counts > 1]),
                " baboons with more than one specimen from\n", length(unique(eth_baboon$PD_ID)), " total baboons)"),
       y=paste0("Worm specimens\n(n=", sum(eth_bab_counts[eth_bab_counts > 1]), " total)"),
       fill="Barcode") +
  guides(x =  guide_axis(angle = 45))
SavePlots(baboon_barcodes_p , plot_path = output_dir, name = "ETH_BaboonBarcodes.png", height = 4, width=6.5)

ClusterPlots(df = dplyr::filter(metadata, country == "Ethiopia" & host == "Baboon (Papio Anubis)"),
            pairwise_df = SubsampleMerge(dplyr::filter(metadata, country == "Ethiopia" & host == "Baboon (Papio Anubis)"), split_variable = NULL),
            output_name =  "ETH_BaboonRelatednessMatrix")

# Work in progress - DAPC analyses
# vcf <- vcfR::read.vcfR(paste0(output_dir, "/", batch_name, "_jointHaploidFilterWithSingletons.vcf.gz"))
# eth_samples <- dplyr::filter(metadata, country == "Ethiopia")
# eth_vcf <- vcf[sample=eth_samples$sample]
# eth_samples <- eth_samples[match(colnames(eth_vcf@gt)[-1], eth_samples$sample),]
# my_genind <- vcfR::vcfR2genind(eth_vcf)
# my_genind@pop <- as.factor(eth_samples$host)
# 
# grp <- find.clusters(my_genind, max.n.clust=40)
# table(pop(my_genind), grp$grp)
# dapc1 <- dapc(my_genind, pop=pop(my_genind), n.pca=20, n.da=6)
# ggsave(filename=paste(output_dir, "CAM-CHD_DAPCByPop.png", sep="/"),
#        plot = call(scatter(dapc1, scree.da=FALSE,  legend=TRUE, col=c("red", "blue"), solid=.4)), 
#        device = "png")

################################################################################
# Mali
################################################################################
mali_sub  <- dplyr::filter(metadata, country == "Mali" & year > 2022)
diff_mali = SubsampleMerge(mali_sub, split_variable = NULL)

ClusterPlots(mali_sub, 
             pairwise_df = diff_mali, 
             output_name="Mali2022+")


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
# foci_order <- c("Northern Jonglei Focal Area", "West Of Nile Focal Area", "Central Focal Area", "East Of Nile Focal Area", "Not provided")
# BarcodeCountPlots(dplyr::filter(metadata, country == "South Sudan") %>%
#                     dplyr::mutate(epi_foci = ifelse(is.na(epi_foci), "Not provided", epi_foci),
#                                   epi_foci = factor(epi_foci, levels=c(names(foci_order), "Not provided"))),
#                   facet_variable = "epi_foci")

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
################################################################################
# This will require updates as the genomics working group comes up with a more detailed analytical plan 

### Subset Chad samples for GA Tech 
chad_counts <- bind_rows(
  GroupByCounts(metadata) %>% mutate(type="Clustered"),
  GroupByCounts(metadata, amplicon_group = "amplicon_barcode_conservative") %>% mutate(type="Conservative"))
SaveTabDelim(chad_counts, paste(output_dir, paste0("Chad_BarcodeCountsYearly-", batch_name), sep="/"))

counts_p <- chad_counts %>% 
  dplyr::mutate(bin = cut(n_specimens, 
                          breaks=c(0,1,5,10,20,50,100,200,300,400,500,1000), 
                          labels=c("1","2-5","6-10","11-20","21-50","51-100","101-200", "201-300", "301-400", "401-500", "> 500"))) %>%
  dplyr::filter(analysis_inclusion == "Included") %>%
  ggplot(aes(x=bin, y=after_stat(count), fill=type)) +
  geom_bar(position="dodge", color="black") +
  labs(y="Unique barcodes", x="Specimens per barcode", fill = "Barcode groups") 
SavePlots(counts_p, output_dir, "CHD_BarcodeBinCountsAll.png", height=4, width=8)
SavePlots(counts_p + facet_grid(year~., scales = "free_y"), output_dir, "CHD_BarcodeBinCountsPerYear.png", height=12, width=8)
rm(chad_counts)  

# This will require updates as the genomics working group comes up with a more detailed analytical plan 
SavePlots(dist_by_sim_plot(dplyr::filter(metadata, country == "Chad"),
                                         country_diff['Chad'][[1]]), 
          output_dir, "CHD_simDist.png", height=5, width=10)

relatedness_summary <- PairwiseWeightedSummary(country_diff["Chad"][[1]], by_year = T)
SavePlots(PairwiseProportionHeatMap(relatedness_summary$weighted_diff),
          output_dir, "CHD_pairwiseHeatmap.png", height=6, width=10)

# Singleton distribution
spatial_singeltons <- dplyr::filter(metadata, country == "Chad" & frequency == 1 & analysis_inclusion == "Included", year >2017) %>%
  ggplot(aes(x=as.numeric(gps_e), y=as.numeric(gps_n), color = confidence, shape = as.character(year))) +
  geom_point() +
  labs(x="Latitude", y="Longitude", shape = "Year", color="Proportion missing")
SavePlots(spatial_singeltons,
          output_dir, "CHD_singletonSpatial.png", height=4, width=6)


# Admin specific comparison
metadata <- tidyr::unite(metadata, "admin_key", admina:admind, remove = F, sep = "-")
metadata <- dplyr::mutate(metadata, admina = gsub("\\ .*", "",admina))
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
 

# Working in progress - Determine ranges and villages that may be likely to have an outbreak
# chad_gps_pairs <- unique(dplyr::select(country_diff['Chad'][[1]], gps_id_pair, meters))
# 
# barcode_range_summary <- as.data.table(dplyr::filter(country_diff['Chad'][[1]], !grepl("_", amplicon_barcode_pair) & !is.na(meters)) %>%
#                                          dplyr::group_by(amplicon_barcode_pair, year_pair) %>%
#                                          summarise(n_pairs = n(),
#                                                    mean_dist = mean(meters), 
#                                                    median_dist = median(meters), 
#                                                    sd_dist=sd(meters)))


# Working in progress - associating changes of new barcodes to case changes
# case_data <- data.table::fread("/mnt/data/guinea_worm/metadata/chad_vas_2020_2024.txt") 
# names(case_data) <- c("year", "admina", "adminb", "adminc", "admind", "surveillance_level", "village_1_plus", "human_case", "animal_infection", "abate", "dog_tether")
# admind_change <- data.table::fread("/mnt/data/guinea_worm/metadata/chad_change_trends_by_snu4_by_year.csv")[-1] %>%
#   tidyr::unite("admin_key", starts_with("snu"), remove = F, sep = "-") 
# names(admind_change) <- c("admin_key", "admina", "adminb", "adminc", "admind", "year", "n_infections", "infection_category", "delta", "delta_category") 
# 
# 
# specimens_wide <- dplyr::filter(metadata, country =="Chad" & !is.na(admina)) %>%
#   group_by(year, admin_key) %>%
#   count()  %>%
#   tidyr::pivot_wider(names_from = year, values_from = n) 
# 
# infection_villages <- dplyr::filter(case_data, animal_infection > 0 | human_case > 0) %>%
#   dplyr::mutate(total_cases = animal_infection + human_case) %>%
#   dplyr::select(year, total_cases, starts_with("admin")) %>%
#   dplyr::arrange(year) %>%
#   tidyr::pivot_wider(names_from = year, values_from = total_cases) %>%
#   dplyr::mutate(infected_years = rowSums(!is.na(select(., starts_with("20"))))) %>%
#   tidyr::unite("admin_key", admina:admind, remove = F, sep = "-")
# point_villages <- dplyr::filter(infection_villages, infected_years == 1) %>%
#   tidyr::pivot_longer(cols = starts_with("20"), names_to = "year") %>%
#   dplyr::filter(!is.na(value)) %>%
#   dplyr::filter(admin_key %in% unique(metadata$admin_key))
# point_vector <- setNames(point_villages$year, point_villages$admin_key)
# 
# columns <- c("year", "amplicon_barcode", "admin_key", "admina", "adminb", "adminc", "admind")
# singletons <- dplyr::filter(metadata, frequency == 1) %>% dplyr::pull(amplicon_barcode)
# 
# 
# admin_barcode_emergence <- dplyr::bind_rows(lapply(specimens_wide$admin_key, FirstBarcodeObservation))
# emergence_ranges <- dplyr::bind_rows(lapply(specimens_wide$admin_key, NearbyBarcodeRanges))
# 
# 
# admin_summary <- admin_barcode_emergence %>%
#   dplyr::group_by(admin_key, barcode_emergence, year) %>%
#   summarize(n_samples = sum(n_specimens)) %>% ungroup() %>% 
#   group_by(admin_key, year) %>%
#   mutate(n_specimens = sum(n_samples),
#          group_proportion = n_samples/n_specimens) %>%
#   dplyr::select(-n_samples) %>%
#   tidyr::pivot_wider(names_from = barcode_emergence, values_from = group_proportion) %>%
#   dplyr::left_join(admind_change, .)



