#######################################################################
# Check Guinea worm metadata
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: 2021/11 
#######################################################################

#######################################################################
# set-up
#######################################################################
for(p in c('data.table', 'dplyr', 'tidyr')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p, repos =  "https://cloud.r-project.org", dependencies = T )
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# set global options
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE)


#######################################################################
# mapping dictionaries
#######################################################################
# set reusable dataframes
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
  c("Unknown", "UNK"),
  c("Unknown", "Unknown"))
names(host_codes) <- c("host", "host_code")

#######################################################################
# functions
#######################################################################
metadata_minor_reformat <- function(df){
  names(df) <- gsub(" ", "_", tolower(names(df)))  
  # provide relevant columns for matching
  df <- df %>%
        dplyr::mutate(vassar_worm = gsub(".*_", "", genomics_sample_id),
                      host = gsub(",domestic", ", domestic", host))
  return(df)
}   


# steps from https://stackoverflow.com/questions/30879429/how-can-i-convert-degree-minute-sec-to-decimal-in-r  
check_gps_degrees <- function(df){
  # check samples for characters commonly associated with degree GPS coordinates
  degrees <- dplyr::filter(df, grepl("\'", gps_n)) %>%
    dplyr::rename(degree_n = gps_n, degree_e = gps_e)
  
  if(nrow(degrees) > 0){
    print(paste("Found", nrow(degrees), "samples that require GPS reformatting."))
    degrees_tmp <- dplyr::select(degrees, degree_n, degree_e) %>% unique() %>%
      tidyr::separate(degree_n, c("N_D", "N_Min", "N_SInt","N_SDec"), remove=F) %>%
      tidyr::unite("N_Sec", "N_SInt","N_SDec", sep = ".") %>%
      tidyr::separate(degree_e, c("E_D", "E_Min", "E_SInt","E_SDec"), remove=F) %>%
      tidyr::unite("E_Sec", "E_SInt","E_SDec", sep=".") %>%
      mutate(across(!starts_with("degree"),
                    ~ as.numeric(as.character(.)))) %>%
      mutate(gps_n = as.character(N_D + N_Min/60 + N_Sec/60^2),
             gps_e = as.character(E_D + E_Min/60 + E_Sec/60^2)) %>%
      dplyr::select(degree_n, gps_n, degree_e, gps_e) %>%
      dplyr::inner_join(., degrees)
    
    # check merging results in non-duplicated samples
    original_length <- nrow(df)
    df <- dplyr::bind_rows(
      dplyr::filter(df, !genomics_sample_id %in% degrees$genomics_sample_id),
      degrees_tmp)
    if(nrow(df) != original_length){
      exit("Effort to convert degrees resulted in larger metadata file.\n
      Check metadata file and retry.")
    } 
   } else {
    print("No entries in the metadata file contain GPS coordinates as degrees.")
    df <- df
  }  
  return(df)  
}

check_sample_metadata <- function(df, vcf_samples = vcf_samples){
  # identify worm numbers that are in the VCF file but not the provided metadata
  missing_meta <- dplyr::filter(vcf_samples, !vassar_worm %in% df$vassar_worm)
  
  if (nrow(missing_meta) > 0){
    print(paste("Found", nrow(missing_meta), "sequences sample that does not match the metadata file.\n
              Parsing the sequencing name to fill high level metadata."))
    name_suffix <- gsub("_.*", "", missing_meta$sample)
    fill_missing <- bind_cols(
      genomics_sample_id = missing_meta$sample,
      vassar_worm = missing_meta$vassar_worm,
      country_code = substr(name_suffix, 1, 3),
      host_code = sapply(name_suffix, function(i) ifelse(nchar(i) < 10, "Unknown", substr(i, 4, 6))),
      year = as.numeric(gsub("[[:alpha:]]", "", name_suffix)),
    )
    fill_missing <- dplyr::left_join(fill_missing, country_codes) %>%
      dplyr::left_join(., host_codes) %>%
      dplyr::select(-host_code, -country_code)
    
    complete_meta <- dplyr::filter(df, vassar_worm %in% df$vassar_worm)
    df <- dplyr::bind_rows(complete_meta, fill_missing) 
  } else{
    print("All worm numbers in the VCF file are represented in the the metadata file.")
    df <- df
  }
  names(df) <- gsub(" ", "_", tolower(names(df)))
  return(df)
}

# add column for which samples have been excluded due to sample coverage
merge_sequencing_quality <- function(df, excluded_samples){
  excluded_worm_numbers <- excluded_samples %>%
    dplyr::mutate(vassar_worm = gsub("^[^_]*_|.batch.*", "", lab_id)) 
  df <- dplyr::left_join(df, excluded_worm_numbers)
  return(df)
}

run_metadata_checks <- function(df, vcf_samples, excluded_samples){
  df <- metadata_minor_reformat(df) %>%
    check_sample_metadata(., vcf_samples) %>%
    check_gps_degrees(.)

  #few final format fixes
  df <- df %>%
    dplyr::mutate(gps_n = gsub(" ", "", gsub(",", ".", gps_n)),
                  gps_e = gsub(" ", "", gsub(",", ".", gps_e)))
  return(df)
}





