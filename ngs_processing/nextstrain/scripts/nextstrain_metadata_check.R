#######################################################################
# Assigns colors for genetic clusters for the GW NextStrain instance. 
# Author: Jessica Ribado, Institute for Disease Modeling
# Date: 2021/11
# Updated: 2024/04 
#######################################################################

#######################################################################
# set-up
#######################################################################
for(p in c('data.table', 'dplyr', 'tidyr', 
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
complete_metadata <- snakemake@input[['barcode_metadata']]
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
    c('Cameroon', 7.3697, 12.3547),
    c("Central African Republic", 6.6194, 20.9367)) 
names(country_gps) <- c("case_gps", "gps_n", "gps_e")


################################################################################
# set colors

# countries
country_colors <- c('#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377',
                    shades::saturation(c('#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377'), shades::scalefac(0.60)))
countries <- c("Chad", "Ethiopia", "Mali", "South Sudan", "Cameroon", "Angola", "Niger", "Sudan", "Cote d'Ivoire", "Central African Republic", "Burkina Faso", "Ghana")
country_colors <- setNames(country_colors, countries)


# genetic groups
gene_base <- data.frame(
    category = "original_barcode",
    value = seq(1,11),
    colors = c("#FF7F00", "#CD0000", "#008B8B", "#7A378B", "#528B8B",
     "#FF6347", "#8DEEEE", "#EEA2AD", "#B4EEB4", "#2F4F4F", "#CDAA7D")
)

amplicon_base <- c("#FF3200FF",  "#D12600", "#DB6A00", "#F9AB0EFF", "#FFED00FF", 
                   "#B2FF2E", "#00AD00", "#019875FF", "#1BB6AFFF", "#32B2DAFF",
                   "#0076BBFF", "#005B94", "#1E2085", "#610052", "#953272", 
                   "#C70E7BFF", "#FC6882FF", "#FF847CFF")
kinship_base <- c("#C62828", "#F44336", "#9C27B0", "#673AB7", "#3F51B5", 
    "#2196F3", "#006064", "#009688", "#4CAF50", "#8BC34A", "#FFEB3B", 
    "#FF9800", "#FFFFFF")
progeny_base <- c("#CC3D24", "#F3C558", "#6DAE90", "#30B4CC", "#004F7A")
reduced_base <- setNames(c("#D9C6B8", "#C2B0A3", "#836F65", "#52271CFF"), 
                         c("Observed once", "Observed in < 5 samples", "Observed in < 10 samples", "Observed in < 20 samples")) 
unaffil_base <- c("#96725B", "#7D4C29", "#613922", "#000000")

amplicon_colors <- setNames(c(amplicon_base, 
                                shades::saturation(amplicon_base, shades::scalefac(0.50)),
                                shades::saturation(amplicon_base, shades::scalefac(0.20))),                                
                                seq(1, length(amplicon_base)*3))
amplicon_colors <- c(amplicon_colors, reduced_base)

# relevant columns
ideal_columns <- c("sample", "host", "year", "sampledate", "country", 
                   "case_gps", "gps_e", "gps_n", "sample_date",
                   "original_barcode", "amplicon_barcode", "amplicon", 
                   "kinship_group", "kinship", "progeny_group", "progeny", 
                   "original_protocol", "microsatellite_protocol")


#######################################################################
# functions
#######################################################################
check_missing_gps <- function(meta_simplified){
  meta_simplified <- dplyr::mutate(meta_simplified,
    gps_n = ifelse(gps_n == "", NA, gps_n),
    gps_e = ifelse(gps_e == "", NA, gps_e))
  missing_gps <-dplyr::filter(meta_simplified, is.na(gps_n) | is.na(gps_e)) %>%
    dplyr::select(-gps_e, -gps_n) %>%
    dplyr::left_join(., dplyr::mutate(country_gps, country = case_gps)) %>%
    dplyr::mutate(gps_e = as.numeric(gps_e), gps_n = as.numeric(gps_n))
  return(missing_gps)
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
        tibble::rowid_to_column("case_gps") %>%
        dplyr::mutate(case_gps = as.character(case_gps))
    metadata <- dplyr::left_join(metadata, gps_coords)
    return(metadata)    
}


geo_tsv <- function(gps_coords){
   total_gps <- rbind.data.frame(
       cbind(Tag = "country", country_gps),
       cbind(Tag = "case_gps", rbind(country_gps, gps_coords))
    )
}



#######################################################################
# run
#######################################################################

#######################################################################
# read in metadata
print(complete_metadata)
metadata <- data.table::fread(complete_metadata, sep="\t", quote="") %>%
  dplyr::filter(excluded_for_analysis == "No") 

filt_columns <- names(metadata)[names(metadata) %in% ideal_columns]
meta_simplified <- dplyr::select(metadata, all_of(filt_columns))

# fill in missing GPS coordinates with country centroid coordinates 
missing_gps_fixed <- check_missing_gps(meta_simplified) 
original_gps <- geo_tag(dplyr::filter(meta_simplified, !sample %in% missing_gps_fixed$sample))
meta_simplified <- dplyr::bind_rows(original_gps, missing_gps_fixed)
if(nrow(meta_simplified) != nrow(metadata)){
  stop("Checking GPS for all successfully sequenced samples produced additional samples. Manually debug issue.")
}

#######################################################################
# assign variable colors

# assign amplicon coloring
ns_colors <- data.table(category ="amplicon", value = names(amplicon_colors), colors = amplicon_colors) 
if(length(unique(meta_simplified$amplicon)) >= length(amplicon_colors)-1){
  print("Amplicon barcodes are within the default color palette. Low frquency observed barcodes (> 10 specimens) will be displayed as individual colors.")
} else{
  print("There are too many colors for the default color palette. Will condense amplicons observed in less than 20 specimens to a single color.")
  freq_11to20 <- unique(dplyr::filter(metadata, frequency > 11 & frequency < 20) %>% .[['amplicon']])
  meta_simplified <- dplyr::mutate(meta_simplified, amplicon = ifelse(amplicon %in% freq_11to20, "Observed in < 20 samples", as.character(amplicon)))
  if(length(unique(meta_simplified$amplicon)) < length(amplicon_colors)){
    print("Collapsing barcodes observed in less than 20 specimens meets default palette length. ")
  } else{
    print("Collapsing barcodes observed in less than 20 specimens still exceeds default palette length.
          High frequency barcode will be colored while lower frequency barcodes will default to grey.
          Improvements to barcode color assignment must be improved manually either in the code or the output color file.")
    ns_colors <- rbind.data.frame(ns_colors,
      cbind(type="amplicon", 
            amplicon = unique(meta_simplified$amplicon)[!unique(meta_simplified$amplicon) %in% color_df$amplicon], 
            color=NA)
      ) 
  }
}

# check for original barcodes 
if("original_barcode" %in% names(meta_simplified)){
    print("Barcode sets with the original protocol provided. ")
} else {
   original_barcode <- read.delim(snakemake@params[['old_barcodes']], sep="\t", header = T)
   meta_simplified <- dplyr::left_join(meta_simplified, original_barcode)
}
meta_simplified <- dplyr::mutate(meta_simplified, 
  original_barcode = as.numeric(meta_simplified$original_barcode),
  original_protocol = ifelse(!is.na(original_barcode), "Sequenced", "Not sequenced"))
ns_colors <- rbind(ns_colors, dplyr::filter(gene_base, value %in% meta_simplified$original_barcode))


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


#######################################################################
# format and save files
meta_simplified <- dplyr::mutate(meta_simplified,
  strain = sample,
  year = as.numeric(year),
  sample_date = ifelse(is.na(sample_date) | sample_date == "", paste0("01/01/", year), sample_date), 
  sample_date = as.Date(sample_date, format = "%m/%d/%Y")) %>%
  dplyr::rename(name = sample, date = sample_date) %>% unique()

# save files
write.table(dplyr::select(meta_simplified, -gps_e, -gps_n), file = out_meta, sep="\t", quote=F, row.names=F) 
write.table(ns_colors, file = out_cols, sep="\t", quote=F, col.names = F, row.names=F) 
write.table(geo_tsv(unique(dplyr::select(meta_simplified, case_gps, gps_n, gps_e))), 
                    file = out_geo, sep="\t", quote=F, col.names = F, row.names=F)



