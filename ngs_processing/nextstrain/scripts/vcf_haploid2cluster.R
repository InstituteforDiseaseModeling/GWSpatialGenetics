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
            'ggsci', 'ggthemes', 'jcolors', 'shades')){
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
metadata <- snakemake@input[['metadata']]
hap_vcf  <- snakemake@input[['haploid_vcf']]
out_meta <- snakemake@output[['meta_clust']]
out_cols <- snakemake@output[['meta_color']]
out_geo  <- snakemake@output[['geo_info']]

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
    c('Cameroon', 12.3547, 7.3697)) 
names(country_gps) <- c("case_gps", "gps_n", "gps_e")

# set colors 
original_bc_colors <- data.frame(
    category = "original_barcode",
    value = seq(1,11),
    #colors = c(ggsci::pal_futurama("planetexpress")(11)[-c(9)], "#7E6148B2")
    colors = c("#FF7F00", "#CD0000", "#008B8B", "#7A378B", "#528B8B",
     "#FF6347", "#8DEEEE", "#EEA2AD", "#B4EEB4", "#2F4F4F", "#CDAA7D")
)

cluster_base <- jcolors::jcolors("rainbow")

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
        id = group %>% grpid,
        identical = ifelse(frequency == 1, "Observed once",
            ifelse(frequency < 5, "Observed in < 5 samples",
                ifelse(frequency < 10, "Observed in < 10 samples", id))))

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

clust_colors <- function(meta_simplified){
    clust_categories <- unique(meta_simplified$identical)
    clust_max <- max(as.numeric(clust_categories), na.rm=T)

    if(clust_max <= length(cluster_base)){
        colors <- cluster_base
    } else if(clust_max <= length(cluster_base)*2){
        colors <- c(cluster_base, as.vector(saturation(cluster_base, 0.5)))
    } else{
        exit(paste0("Too many colors (>",  length(cluster_base)*2, ") for automatic coloring, requires manual color palette updates. Exiting."))
    }

    cluster_colors <- data.frame(
        category = "identical",
        value = c(seq(1, clust_max), "Observed once", "Observed in < 5 samples", "Observed in < 10 samples"),
        colors = c(colors[1:clust_max], "#999999", "#666666", "#333333")
    ) 
    print(cluster_colors %>% head())
    return(cluster_colors)
} 

#######################################################################
# run
#######################################################################
# identify unique clusters
vcf <- vcfR::read.vcfR(hap_vcf, verbose = FALSE)
vcf_clust <- gt2barcode(vcf)

# add to metadata
print(metadata)
metadata <- read.delim(metadata, sep="\t") 
names(metadata) <- tolower(names(metadata))

if("original_barcode" %in% names(metadata)){
    print("Barcode sets with the original protocol provided. ")
    metadata$original_barcode <- as.numeric(metadata$original_barcode)
} else{
    print("Barcode sets with the original protocol were not in the metadata file provided. Setting to NA.")
    metadata$original_barcode <- NA
}

meta_cluster <- dplyr::inner_join(metadata, vcf_clust) %>%
    left_join(., geo_tag(metadata)) 
write.table(meta_cluster, file = gsub("_clusters", "_allInfo", out_meta), sep="\t", quote=F, row.names=F) 

meta_simplified <- dplyr::select(meta_cluster, sample, host, year, sampledate, country, case_gps, original_barcode, id, identical) %>%
    dplyr::mutate(
        strain = sample,
        case_gps = ifelse(is.na(case_gps), country, case_gps), 
        sampledate = ifelse(sampledate == "", paste0("01/01/",year), sampledate), 
        sampledate = as.Date(sampledate, format = "%m/%d/%Y")) %>%
    dplyr::rename(name = sample, date = sampledate) %>% unique()
write.table(meta_simplified, file = out_meta, sep="\t", quote=F, row.names=F) 

nextidentical_cols <- rbind.data.frame(
    original_bc_colors, 
    clust_colors(meta_simplified)
) 
write.table(nextidentical_cols, file = out_cols, sep="\t", row.names=F, col.names=F, quote=F)