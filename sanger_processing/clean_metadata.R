# Cleaning metadata
# Jessica Ribado
# March 2020
# 
# Note: 2015 data does not separate quarters from villages, and there's not easy way to do this programatically since these columns are indicated by coloring on excel. This shouldn't change anything other than adding columns that won't merge. 
################################################################################
# load libraries
################################################################################ 
require("readxl")
require(tidyverse)
 
################################################################################
# load and format national line list data 
# will provide GPS coordinates for any worms, be used to check any inconsistencies between years for the same region
################################################################################
natLine_dir <- '/home/jribado/Dropbox (IDM)/GW modeling - data sharing repository/Surveillance_data_national_line_lists'
# leave out 2018 since the file is structured differneyl and requires different processing
natLine_files <- list(natLine_2015='National Line List_Dec2015_01.15.16.xlsx', 
                      natLine_2016='National Line List - Decembre 2016.xlsx',
                      natLine_2017='1 - National Line List - Dec 2017 - Master 13.01.18.xlsx') 
natLine_list <- lapply(setNames(names(natLine_files), names(natLine_files)), function(i){
  natLine_tmp <- readxl::read_excel(paste(natLine_dir, natLine_files[[i]], sep="/"), sheet=1, col_names = TRUE)
  names(natLine_tmp) <- natLine_tmp[1,]
  natLine_tmp <- natLine_tmp[3:nrow(natLine_tmp), c("District", "Nom de Villages", "Latitude", "Longitude")]
  natLine_tmp <- dplyr::rename(natLine_tmp, "Village" = `Nom de Villages`) %>%
    dplyr::mutate(year = gsub("natLine_", "", i), 
                  Village = gsub(" II", " 2", Village))
  assign(i, natLine_tmp)
})
# edit 2018 data separately
natLine_18 <- readxl::read_excel(paste(natLine_dir, 'National Line List _ Dec2018__final.xlsx', sep="/"), sheet=1, col_names = TRUE)
names(natLine_18) <- natLine_18[2,]
natLine_18 <- natLine_18[3:2754,c("District", "Village", "Quartier", "Latitude", "Longitude")]
natLine_18 <- dplyr::filter(natLine_18, is.na(Quartier)) %>% dplyr::select(-Quartier)
natLine_18$year <- "2018"
# merge together
natLine_merge <- rbind(do.call(rbind, natLine_list), natLine_18) %>%
  dplyr::mutate(Latitude = round(as.numeric(Latitude), 5), Longitude = round(as.numeric(Longitude), 5)) %>%
  dplyr::arrange(Latitude, Longitude) %>% unique() %>%
  tidyr::spread(year, year)
names(natLine_merge) <- tolower(names(natLine_merge)) 
# write table
# write.table(natLine_merge, '/home/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Guinea Worm Genetics in Chad/sample_identifiers/20203_natLine_merging.txt', sep="\t", row.names = F, quote = F) 


################################################################################
# load and format case data
################################################################################
case_dir <- '/home/jribado/Dropbox (IDM)/GW modeling - data sharing repository/'
human_casesList <- lapply(seq(1,4), function(i){
  readxl::read_excel(paste(case_dir, 'Surveillance_data_humans/HUMANS_Line_Listing of Cases_Updated_18oct182018_deidentified.xlsx', sep="/"), sheet=i, col_names = TRUE)
})
dog_cases15to18 <- readxl::read_excel(paste(case_dir, 'Surveillance_data_worms/worms_15_18.xlsx', sep="/"), sheet=1, col_names = TRUE)
dog_cases19 <- read.delim(paste(case_dir, 'Surveillance_data_worms/VG Chez Les Animaux_2020_14_1_20_deidentified.csv', sep="/"), sep="\t")
names(dog_cases19) <- dog_cases19[1,]

# subset relavent columns and rename columns for merging
case_colNames <- c("worm_number", "district", "village", "latitude", "longitude", "year", "month", "day")
french_months <- cbind.data.frame(
  french_month = c("janv", "oct", "août", "févr", "mars", "mai",  "avr",  "nov",  "sept", "juin", "juil", "déc", "Dec", "Nov", "Aug", "Oct"),
  month = c(1, 10, 8, 2, 3, 5, 4, 11, 9, 6, 7, 12, 12, 11, 8, 10)
)

# humans 
human_sub <- dplyr::bind_rows(lapply(seq(1,4), function(i){
  if(i %in% c("1", "2")){
    dplyr::select(human_casesList[[i]], Confirmed, District, `Village du Detection`, `gps latitude`, `gps longitude`, `Date d'Emergence`) %>% 
      tidyr::separate(`Date d'Emergence`, c("year", "month", "day"))
  } else{
    year <- ifelse(i == "3", "2016", "2015")
    dplyr::select(human_casesList[[i]], Confirmed, District, `Village du Detection`, `gps latitude`, `gps longitude`) %>% 
     dplyr::mutate(year = year, month ="NA", day="NA")
  }
}))
names(human_sub) <- case_colNames

# dogs
dog15to18_sub <- dplyr::select(dog_cases15to18, worm_number, district_name, village_name, `GPS_N`, `GPS_E`, `date_emerge_sas`)  %>% 
  tidyr::separate(`date_emerge_sas`, c("year", "month", "day"))
names(dog15to18_sub) <- case_colNames
# removes duplicate column with same name that gives an error
dog19_sub <- dplyr::select(dog_cases19[-1,-31], `N. de ver `, District, Village, `Date d'émergence`)%>% 
  tidyr::separate(`Date d'émergence`, c("day", "french_month", "year")) %>%
  dplyr::left_join(., french_months) %>%
  dplyr::mutate(year = "2019", latitude = NA, longitude = NA)
dog19_sub <- dplyr::select(dog19_sub, `N. de ver `, District, Village, latitude, longitude, year, month, day)
names(dog19_sub) <- case_colNames
# merge datasets for all cases available
cases_merged <- rbind(human_sub, dog19_sub, dplyr::filter(dog15to18_sub, is.na(latitude) | !grepl("[A-Z]", latitude))) %>%
  dplyr::mutate(latitude = round(as.numeric(latitude), 5), 
                longitude = round(as.numeric(longitude), 5)) %>% unique()
# write.table(cases_merged, '/home/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Guinea Worm Genetics in Chad/sample_identifiers/20203_cases_merged.txt', sep="\t", row.names = F, quote = F)


################################################################################
# Combine national line and case data
################################################################################
# for worms with gps coorindates, merge with the national line to confirm which ones are correct

merged_withGPS <- merge(dplyr::filter(cases_merged, !is.na(latitude) & !is.na(longitude)),
                        natLine_merge, by=c("district", "village", "latitude", "longitude"), all.x=T)
merged_missingGPS <- merge(dplyr::filter(cases_merged, is.na(latitude) & is.na(longitude)) %>% dplyr::select(-latitude, -longitude), 
                           natLine_merge, by=c("district", "village"), all.x=T)
merged_meta <- rbind(dplyr::select(merged_withGPS, c(case_colNames, `2015`, `2016`, `2017`, `2018`)),
                     dplyr::select(merged_missingGPS, c(case_colNames, `2015`, `2016`, `2017`, `2018`))) %>%
  dplyr::add_count(worm_number, year)
# write.table(merged_meta, '/home/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Guinea Worm Genetics in Chad/sample_identifiers/202004_natLineCasesDuplicates.txt', sep="\t", row.names = F, quote = F)
