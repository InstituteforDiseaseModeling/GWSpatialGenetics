# guinea worm functions - running genetic pairwise similarity
# not needed as of 05/2020 - these functions are for protein not nucleotide similarity

################################################################################
# genetic similarity functions
################################################################################
df2barcodeList <- function(df){
  df[df == "-"] <- "N"
  barcodes <- as.list(df[,2])
  names(barcodes) <- df[,1]
  return(barcodes)
}

simMat <- function(barcode_list, cores, names_list){
  sim_mat <- parSeqSim(barcode_list, cores = cores, type = "local")
  colnames(sim_mat) <- names_list
  rownames(sim_mat) <- names_list
  return(sim_mat)
}

simMat2df <- function(simMat){
  # keep only the lower part of the matrix to keep only one correlation per pair
  simMat[upper.tri(simMat)] <- NA
  pw_df <- reshape2:::melt.matrix(simMat, varnames = c('worm1', 'worm2'), na.rm = TRUE)
  # filter out self-correlations
  pw_df <- dplyr::filter(pw_df, worm1 != worm2) %>%
    mutate_if(is.factor, as.character)
  return(pw_df)
}

simRun <- function(barcode_df, threads){
  tmp_bc  <- df2barcodeList(barcode_df)
  tmp_mat <- simMat(tmp_bc, 16 , names(singletonInc_bc))
  sim_df  <- simMat2df(tmp_mat)
  return(sim_df)
}

################################################################################
# new metadata columns to use in epi-spatial models
################################################################################
pairwise_comp <- function(sim_df, metadata_df, barcode_df){
  bc_to_md <- merge(metadata_df, barcode_df, by="worm")
  sim_merge <- dplyr::left_join(bc_df, dplyr::rename(bc_to_md, worm1 = worm), by="worm1")
  sim_merge <- dplyr::left_join(sim_merge, dplyr::rename(bc_to_md, worm2 = worm), by="worm2")
  sim_merge <- dplyr::mutate(sim_merge,
    bc_match = ifelse(full_barcode.x == full_barcode.y, "True", "False"),
    yr_match = ifelse(year.x == year.y, "True", "False"),
    host_match = ifelse(host_number.x == host_number.y, "True", "False")
    #distance = ifelse(!is.na(latitude.x) & !is.na(longitude.x) & !is.na(latitude.y) & !is.na(longitude.y),
    #  geosphere::distm(c(latitude.x, longitude.x), c(latitude.y, longitude.y), fun=distHaversine), "NA")
  )
}
