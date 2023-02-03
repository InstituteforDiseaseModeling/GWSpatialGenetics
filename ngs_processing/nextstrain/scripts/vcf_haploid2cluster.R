#######################################################################
# Creates barcode clusers based on haploid VCF file for spatial analyses
# Note: Will be updated with a new 
# Author: Jessica Ribado, Institute for Disease Modeling
# Date: 2021/11 
#######################################################################

#######################################################################
# set-up
#######################################################################
for(p in c('vcfR', 'data.table', 'dplyr', 'tidyr', 
            'ggsci', 'ggthemes', 'shades', 'paletteer')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p, repos =  "https://cloud.r-project.org", dependencies = T )
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# set global options
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE)

################################################################################
# load variables from snakemake
metadata <- snakemake@input[['metadata']]
hap_vcf  <- snakemake@input[['haploid_vcf']]
out_meta <- snakemake@output[['meta_clust']]
out_cols <- snakemake@output[['meta_color']]
out_geo  <- snakemake@output[['geo_info']]

################################################################################
# set reusable dataframes
country_gps <- rbind.data.frame(
    c("Sudan", 12.8628, 30.2176),
    c("South Sudan", 4.85, 31.6),
    c('Chad', 15.4542, 18.7322),
    c('Niger',	17.6078, 8.0817),
    c('Mali',  17.5707,	-3.9962),
    c('Burkina Faso', 12.2383, -1.5616),
    c("Cote d'Ivoire", 7.54, -5.5471),
    c('Ghana', 7.9465, -1.0232),
    c('Ethiopia',  9.145, 40.4897),
    c('Angola', -11.2027, 17.8739),
    c('Cameroon', 12.3547, 7.3697),
    c("Central African Republic", 6.6194, 20.9367)) 
names(country_gps) <- c("case_gps", "gps_n", "gps_e")

country_codes <- rbind.data.frame(
  c("Sudan", "SUD"),
  c("South Sudan", "SSU"),
  c('Chad', "CHD"),
  c('Niger',	"NGR"),
  c('Mali',  "MAL"),
  c('Burkina Faso', "BFA"),
  c("Cote d'Ivoire", "CDI"),
  c('Ghana', "GHA"),
  c('Ethiopia',  "ETH"),
  c('Angola', "ANG"),
  c('Cameroon', "CAM"),
  c("Central African Republic", "RCA")) 
names(country_codes) <- c("country", "country_code")

host_codes <- rbind.data.frame(
  c("Dog", "DOG"),
  c("Human", "HUM"),
  c("Cat, domestic", "CAT"),
  c("Cat, Panthera pardus", "PPD"),
  c("Cat, wild unknown spp.", "WCT"),
  c("Baboon", "BAB"),
  c("Unknown animal", "Unidentified")
)
names(host_codes) <- c("host", "host_code")

################################################################################
# set colors 
gene_base <- data.frame(
    category = "original_barcode",
    value = seq(1,11),
    colors = c("#FF7F00", "#CD0000", "#008B8B", "#7A378B", "#528B8B",
     "#FF6347", "#8DEEEE", "#EEA2AD", "#B4EEB4", "#2F4F4F", "#CDAA7D")
)

amplicon_base <- c("#D12600", "#DB6A00", "#B2FF2E", "#00AD00", "#9CCADE", 
                  "#005B94", "#1E2085", "#610052", "#953272")
kinship_base <- c("#C62828", "#F44336", "#9C27B0", "#673AB7", "#3F51B5", 
    "#2196F3", "#006064", "#009688", "#4CAF50", "#8BC34A", "#FFEB3B", 
    "#FF9800", "#FFFFFF")
progeny_base <- c("#CC3D24", "#F3C558", "#6DAE90", "#30B4CC", "#004F7A")
reduced_base <- setNames(c("#D9C6B8", "#C2B0A3", "#836F65"), 
                         c("Observed once", "Observed in < 5 samples", "Observed in < 10 samples")) 
unaffil_base <- c("#96725B", "#7D4C29", "#613922", "#000000")


#######################################################################
# functions
#######################################################################
grpid <- function(x) match(x, unique(x))

gt2barcode <- function(vcf){
  # get exact genotypes 
  tmp_gt <- data.frame(extract.gt(vcf, element = "GT", return.alleles = TRUE))
  
  tmp_gt[is.na(tmp_gt)] <- "." 
  sequences <- data.frame(sequence = apply(tmp_gt, 2, paste, collapse="")) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::group_by(sequence) %>% 
    dplyr::mutate(group = group_indices(),
                  frequency = n()) %>%
    dplyr::arrange(desc(frequency), group) %>% ungroup() %>%          
    dplyr::mutate(
        amplicon_barcode = group %>% grpid,
        amplicon = ifelse(frequency == 1, "Observed once",
            ifelse(frequency < 5, "Observed in < 5 samples",
                ifelse(frequency < 10, "Observed in < 10 samples", amplicon_barcode))))

    return(sequences)
}

geo_tag <- function(metadata){

    if(!"sample" %in% colnames(metadata)){
        exit("Metadata file missing 'sample' column.")
    } 

    if(!(c("gps_e") %in% colnames(metadata)) | !(c("gps_e") %in% colnames(metadata))){
        exit("Metadata file missing longitude (GPS_E) and/or latitude (GPS_N) columns.")
    } 

    gps_coords <- dplyr::select(metadata, gps_n, gps_e) %>%
        drop_na() %>% unique() %>%
        tibble::rowid_to_column("case_gps") 

    geo_tsv(gps_coords, out_geo)    

    return(gps_coords)    
}

geo_tsv <- function(gps_coords, output_file){
   total_gps <- rbind.data.frame(
       cbind(Tag = "country", country_gps),
       cbind(Tag = "case_gps", rbind(country_gps, gps_coords))
    )
    write.table(total_gps, file = output_file, quote = F, row.names = F, col.names = F, sep="\t")     
}

nextstrain_colors <- function(df, column, base_colors){
  # double the number of available colors
  expanded_colors <- c(base_colors, 
                       saturation(paste0(base_colors, "FF"), scalefac(0.5)),
                       saturation(paste0(base_colors, "FF"), scalefac(0.25)))
  
  col_entries <- dplyr::pull(df, column)
  ungrouped  <- sort(unique(col_entries[grepl("unaffiliated", col_entries)]))
  grouped <- names(sort(table(col_entries[!grepl("unaffiliated|Observed", col_entries)]), decreasing=T))
  if(length(ungrouped) > 0){
    merged_group <- c(setNames(expanded_colors[1:length(grouped)], grouped), 
                      setNames(unaffil_base[1:length(ungrouped)], ungrouped))
  }else{
    merged_group <- setNames(expanded_colors[1:length(grouped)], grouped)
  }

  colors <- cbind(category = column,
    data.frame(colors = c(merged_group,
    reduced_base[names(reduced_base) %in% unique(col_entries)])) %>%
    tibble::rownames_to_column("value"))
  return(colors)
}
  
  
#######################################################################
# run
#######################################################################
# identify unique clusters
vcf <- vcfR::read.vcfR(hap_vcf, verbose = FALSE)
vcf_clust <- gt2barcode(vcf)

# add to metadata
metadata <- read.delim(metadata, sep="\t") %>%
  dplyr::mutate(vassar_worm = gsub(".*_", "", Genomics.Sample.ID)) %>%
  dplyr::rename("SampleDate"="Sample.date")
names(metadata) <- tolower(names(metadata))
metadata[metadata == "." & !is.na(metadata)] <- NA

# check that metadata exists for all files in the VCF
vcf_samples <- colnames(vcf@gt)[-1]
vcf_numbers <- cbind.data.frame(
  sample = vcf_samples,
  vassar_worm = gsub(".*_", "", vcf_samples))
# merge back with metadata
metadata <- dplyr::left_join(metadata, vcf_numbers)

yes_meta <- dplyr::filter(metadata, vassar_worm %in% vcf_numbers$vassar_worm)
missing_meta <- vcf_samples[!vcf_samples %in% yes_meta$sample]
if (length(missing_meta) > 0){
  name_suffix <- gsub("_.*", "", missing_meta)
  fill_missing <- bind_cols(
    sample = missing_meta,
    country_code = substr(name_suffix, 1, 3),
    host_code = sapply(name_suffix, function(i) ifelse(nchar(i) < 10, "Unidentified", substr(i, 4, 6))),
    year = as.numeric(gsub("[[:alpha:]]", "", name_suffix))
  )
  
  fill_missing <- dplyr::inner_join(fill_missing, country_codes) %>%
    dplyr::inner_join(., host_codes) %>%
    dplyr::select(-host_code, -country_code)
  metadata <- dplyr::bind_rows(fill_missing, yes_meta) %>%
    dplyr::filter(sample %in% vcf_samples)
  #print(head(metadata))
  if(nrow(metadata) != length(vcf_samples)){
    print(paste("Metadata available for", nrow(metadata), "samples."))
    print(paste("VCF files contains", length(vcf_samples), "samples prior to filtering."))
    #stop("The metadata columns do not match the number of samples in the VCF after parsing names. Please check input names.")
  }
}


if("original_barcode" %in% names(metadata)){
    print("Barcode sets with the original protocol provided. ")
    metadata <- dplyr::mutate(metadata, 
      original_barcode = as.numeric(metadata$original_barcode),
      original_protocol = ifelse(!is.na(original_barcode), "Sequenced", "Not sequenced"))
    ns_colors <- dplyr::filter(gene_base, value %in% metadata$original_barcode)
} 

# check the number of samples per group to assign new colors for grouping, if needed
min_samples = 5
if("kinship_group" %in% names(metadata)){
  print("Kinship groups are provided.")
  if (length(unique(na.omit(metadata$kinship_group))) > length(kinship_base)*2){
    print(" More kinships than available colors. Assigning kinship groups found in less than 5 samples to a single color. ")
    metadata <- dplyr::group_by(metadata, kinship_group) %>%
      dplyr::mutate(kinship_frequency = n(),
                    kinship = ifelse(kinship_frequency < min_samples | grepl("unaffliated", kinship_group), paste("Observed in <", min_samples, "samples"), kinship_group))
    if(length(unique(na.omit(metadata$kinship))) >  length(kinship_base)*2){
      exit(paste0("Too many colors (>",  length(kinship_base)*2, ") for automatic coloring, requires manual color palette updates. Exiting."))
    }
  } else {
    metadata <- dplyr::group_by(metadata, kinship_group) %>%
      dplyr::mutate(kinship_frequency = n(),
                    kinship = kinship_group)
  }
  kinship_colors <- nextstrain_colors(metadata, "kinship", kinship_base)
  if(exists("ns_colors")){
    ns_colors <- rbind(ns_colors, kinship_colors)
  } else{
    ns_colors <- kinship_groups
  }
  
  if("progeny_group" %in% names(metadata)){
    if (length(unique(na.omit(metadata$progeny_group))) > length(progeny_base)*2){
      print("More parent offspring pairs than available colors.
            Assigning pairs with less than 5 samples per family to a single color. ")
      metadata <- dplyr::group_by(metadata, progeny_group) %>%
        dplyr::mutate(progeny_frequency = n(),
                      progeny = ifelse(progeny_frequency < min_samples, paste("Observed in <", min_samples, "samples"), progeny_group)) 
    } else{
      metadata <- dplyr::group_by(metadata, progeny_group) %>%
        dplyr::mutate(progeny_frequency = n(),
                      progeny = progeny_group)
    }
    ns_colors <- rbind(ns_colors, nextstrain_colors(metadata, "progeny", progeny_base))
  }  
  metadata <- metadata %>% ungroup()
  metadata$microsatellite_protocol <- ifelse(!is.na(metadata$kinship_group), "Sequenced", "Not sequenced")
} 

meta_cluster <- dplyr::inner_join(metadata, vcf_clust) %>%
    left_join(., geo_tag(metadata)) 
write.table(meta_cluster, file = gsub("_clusters", "_allInfo", out_meta), sep="\t", quote=F, row.names=F) 

ideal_columns <- c("sample", "host", "year", "sampledate", "country", "case_gps",
                  "original_barcode", "amplicon_barcode", "amplicon", 
                  "kinship_group", "kinship", "progeny_group", "progeny", 
                  "original_protocol", "microsatellite_protocol")
filt_columns <- names(meta_cluster)[names(meta_cluster) %in% ideal_columns]

meta_simplified <- dplyr::select(meta_cluster, all_of(filt_columns)) %>%
    dplyr::mutate(
        strain = sample,
        case_gps = ifelse(is.na(case_gps), country, case_gps), 
        year = as.numeric(year),
        sampledate = ifelse(sampledate == "", paste0("01/01/",year), sampledate), 
        sampledate = as.Date(sampledate, format = "%m/%d/%Y")) %>%
    dplyr::rename(name = sample, date = sampledate) %>% unique()
write.table(meta_simplified, file = out_meta, sep="\t", quote=F, row.names=F) 

nextidentical_cols <- nextstrain_colors(meta_simplified, "amplicon", amplicon_base)
# check if other genotype categories need to be included in the output file
if(exists("ns_colors")){
  nextidentical_cols <- rbind.data.frame(nextidentical_cols, ns_colors) 
} 

write.table(nextidentical_cols, file = out_cols, sep="\t", row.names=F, col.names=F, quote=F)
