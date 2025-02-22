################################################################################
# Commonly used functions for Guinea worm genomic analyses
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: May 2023, updated February 2025
################################################################################

################################################################################
# General
################################################################################
SaveTabDelim <- function(r_obj, path){
  write.table(r_obj, file = path, quote = F, row.names = F, sep="\t")
}


mutate_cond <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data %>% filter(condition) %>% mutate(...)
  .data
}


################################################################################
# Pairwise distance functions
################################################################################
HaversineDist <- function(gps_coords){
  d0 <- geodist::geodist(gps_coords %>% 
                           dplyr::rename(longitude = "gps_e", latitude = "gps_n"),
                         measure = "haversine")
  colnames(d0) <- gps_coords$gps_id
  row.names(d0) <- gps_coords$gps_id
  
  dist_df <- reshape2::melt(d0) %>% 
    dplyr::rename(gps.x=Var1, gps.y=Var2, meters=value) %>%
    dplyr::filter(gps.x != gps.y) %>%
    rowwise() %>%
    dplyr::mutate(gps_id_pair = paste(min(gps.x, gps.y), max(gps.x, gps.y), sep="_"))
  
  return(dist_df)  
}


################################################################################
# Pairwise metadata merging functions
################################################################################
PairwiseMetadataGroups <- function(df, dt){
  # Uses vectorization of multiple columns to identify matches for groupings instead of dplyr/merging for speed. 
  # Giving odd errors still on indexing values. Skip for now. 
  compare_columns <- c("amplicon_barcode", "year", "host", "emergence_date", "gps_id")
  lookup_tables <- lapply(setNames(compare_columns, compare_columns), function(i){
    setNames(as.character(df[[i]]), df[["sample"]])
  })
  
  for (lookup_name in names(lookup_tables)) {
    lookup_vector <- lookup_tables[[lookup_name]]
    new_col_name <- paste0(lookup_name, "_pair")  # Column name for the matched category
    
    # dt <- dt[, (paste0(lookup_name, ".x")) := lookup_vector[match(sample.x, names(lookup_vector))]]
    # dt <- dt[, (paste0(lookup_name, ".y")) := lookup_vector[match(sample.y, names(lookup_vector))]]
    
    if(lookup_name %in% c("amplicon_barcode", "gps_id", "year")){
      dt <- dt[, (new_col_name) := data.table::fcase(
        lookup_vector[match(sample.x, names(lookup_vector))] == lookup_vector[match(sample.y, names(lookup_vector))], as.character(lookup_vector[match(sample.x, names(lookup_vector))]),
        is.na(lookup_vector[match(sample.x, names(lookup_vector))]) | is.na(lookup_vector[match(sample.y, names(lookup_vector))]), "Missing data",
        default = paste(pmin(lookup_vector[match(sample.x, names(lookup_vector))], lookup_vector[match(sample.y, names(lookup_vector))]),
                      pmax(lookup_vector[match(sample.x, names(lookup_vector))], lookup_vector[match(sample.y, names(lookup_vector))]), sep="_")
      )]
    } else if(lookup_name == "emergence_date"){
      dt <- dt[, (new_col_name) := abs(as.numeric(
        difftime(lookup_vector[match(sample.x, names(lookup_vector))], lookup_vector[match(sample.y, names(lookup_vector))], 
                 unit="days")))]
      dt <- dt[, generation := data.table::fcase(
        between(emergence_date_pair, 310, 410), "Within 10-14 mo",
        default = "Unlikely or Unknown")]
    } else{
      dt <- dt[, (new_col_name) := data.table::fcase(
        lookup_vector[match(sample.x, names(lookup_vector))] == lookup_vector[match(sample.y, names(lookup_vector))], as.character(lookup_vector[match(sample.x, names(lookup_vector))]),
        is.na(lookup_vector[match(sample.x, names(lookup_vector))]) | is.na(lookup_vector[match(sample.y, names(lookup_vector))]), "Missing data",
        default = "Mixed")]
    }
  } 
  return(dt)
}  


SubsampleMerge <- function(df, diff_df = barcode_pairwise, split_variable = "year"){
  
  UniquePairsDiff <- function(df, unique_pairs){
    # merge in pairwise metadata 
    if(typeof(unique_pairs) %in% c("character", "matrix")){
      unique_pairs <- t(unique_pairs)
    } 
    pairs_dt <- as.data.table(unique_pairs)
    setnames(pairs_dt, c("V1", "V2"), c("sample.x", "sample.y"))
    dt_meta <- PairwiseMetadataGroups(df, pairs_dt)

    tmp_same_bc <- dt_meta[!grepl("_", amplicon_barcode_pair),]
    tmp_same_bc <- tmp_same_bc[, `:=`(all = 0, rmNA = 0, missing = 0)]

    tmp_diff_bc <- dt_meta[grepl("_", amplicon_barcode_pair),]
    tmp_diff_bc <- merge(tmp_diff_bc, diff_df[ ,c("amplicon_barcode_pair", "all", "rmNA", "missing")], 
                         by = "amplicon_barcode_pair", all.x=T, allow.cartesian = T)

    pairs_df <- dplyr::bind_rows(tmp_same_bc, tmp_diff_bc)
    rm(pairs_dt, dt_meta, tmp_same_bc, tmp_diff_bc)
    
    # Calculate Haversine distance for the samples in this set, add to dataset
    gps_dist <- HaversineDist(unique(dplyr::select(df, gps_id, gps_e, gps_n)))
    pairs_df <- merge(pairs_df, unique(dplyr::select(gps_dist, gps_id_pair, meters)),
                      by = "gps_id_pair", all.x=T, allow.cartesian = T)

    return(pairs_df)
  }
  
  # Initialize a vector with elements
  if(is.null(split_variable)){
    print("Running all-by-all pairwise for specified samples.")
    unique_pairs <- combn(df[['sample']], 2)
    print(paste("Total pairwise comparisons:", ncol(unique_pairs)))
    sample_pairwise <- UniquePairsDiff(df, unique_pairs)
  } else{
    # specify groups of pairs to run relatedness
    sample_split <- split(df[["sample"]], df[[split_variable]])
    
    if(split_variable == "year"){
      # run all-by-all specimen pairs within a year
      within_years <- parallel::mclapply(names(sample_split), function(i){
        if(length(sample_split[i][[1]]) > 1){
          unique_pairs <- combn(sample_split[i][[1]], 2)
          UniquePairsDiff(df, unique_pairs)
        }
      })
      # check for prior year in data set, run all-by-all specimen pairs with the previous year
      years_to_check = as.numeric(names(sample_split))
      previous_years <- lapply(years_to_check, function(j){
        if((j-1) %in% years_to_check){
          year_pairs <- expand.grid(sample_split[as.character(j-1)][[1]],
                                    sample_split[as.character(j)][[1]])
          
          names(year_pairs) <- c("V1", "V2")
          UniquePairsDiff(df, year_pairs)
        }
      })
    sample_pairwise <- dplyr::bind_rows(within_years, previous_years) 
    rm(within_years, previous_years)
    }
  }
  return(sample_pairwise)
}
  
 

################################################################################
# Summary statistic functions
################################################################################
# https://stackoverflow.com/questions/2748725/is-there-a-weighted-median-function
weighted.median <- function(x, w, q=.5) {
  n <- length(x)
  i <- order(x)
  w <- cumsum(w[i])
  p <- w[n] * q
  j <- findInterval(p, w)
  Vectorize(function(p,j) if(w[n] <= 0) NA else
    if(j < 1) x[i[1]] else
      if(j == n) x[i[n]] else
        if(w[j] == p) (x[i[j]] + x[i[j+1]]) / 2 else
          x[i[j+1]])(p,j)
}


PairwiseWeightedSummary <- function(pairwise_df, by_year = FALSE){
  pairwise_long <- dplyr::filter(pairwise_df, !is.na(rmNA)) %>%
    dplyr::mutate(year_pair = ifelse(grepl("_", year_pair), "Mixed", year_pair)) %>%
    dplyr::group_by(rmNA, !!sym(ifelse(by_year, "year_pair", ""))) %>%
    dplyr::summarise(pairwise_count = n()) 
  
  pairwise_summary <- pairwise_long %>%
    summarize(pairwise_sum = sum(pairwise_count),
          weighted_mean = weighted.mean(1-rmNA, pairwise_count),
          weighted_median = weighted.median(1-rmNA, pairwise_count),
           # Variance formula: sum((weights * (values - weighted_mean)^2)) / sum(weights)
           weighted_variance = sum(pairwise_count * (1-rmNA - weighted_mean)^2) / sum(pairwise_count),
           weighted_sd = sqrt(weighted_variance))
  return(list(weighted_diff = pairwise_long, 
              summary_diff = pairwise_summary))
}



################################################################################
# Deprecated - All-by-all relatedness functions
################################################################################
# Uses a matrix of all samples - slow
PairwiseDifference <- function(gt_mat){
  diff_rmMissing <- ape::dist.gene(gt_mat, method="percentage", pairwise.deletion = TRUE)
  gt_sub <- gt_mat
  gt_sub[is.na(gt_sub)] <- "X"
  diff_all <- ape::dist.gene(gt_sub, method="percentage")
  diff_df <- dplyr::inner_join(reshape2::melt(LowerMatrix(diff_all), na.rm=T) %>% rename("All"=value),
                               reshape2::melt(LowerMatrix(diff_rmMissing), na.rm=T) %>% rename("rmNA"=value))
  
  # identify the indices of missing positions to count the number of excluded positions in a pair of sample
  na_pos <- lapply(setNames(rownames(gt_mat),rownames(gt_mat)), 
                   function(i) as.vector(which(is.na(gt_mat[i,]))))
  diff_df$missing = MissingCountsV(na_pos, diff_df$Var1, diff_df$Var2)
  diff_df <- dplyr::mutate_if(diff_df, is.numeric, round, digits=3)
  return(diff_df)
} 

LowerMatrix <- function(matrix){
  if(class(matrix) != "matrix"){
    matrix <- as.matrix(matrix)
  }
  matrix[upper.tri(matrix, diag = FALSE)] <- NA
  return(matrix)
}

MissingCounts <- function(pos_list, sample1, sample2){
  pairwise_missing <- unlist(c(pos_list[sample1], pos_list[sample2]))
  return(length(pairwise_missing))
}
MissingCountsV <- Vectorize(MissingCounts, vectorize.args = c("sample1", "sample2"))


AllPairwise <- function(vcf){
  # Extract genotypes from the VCF file
  tmp_gt <- data.frame(extract.gt(vcf, element = "GT", return.alleles = TRUE))
  gt <- data.frame(t(tmp_gt))
  names(gt) <- row.names(tmp_gt)
  row.names(gt) <- names(tmp_gt)
  genetic_diff <- PairwiseDifference(gt)
  
  diff_df <- unique(dplyr::left_join(genetic_diff, pairwise_dist)) %>%
    mutate_cond(Var1==Var2, All=NA, rmNA=NA, missing=NA, meters=NA) %>%
    dplyr::rename("sample.x"=Var1, "sample.y"=Var2)
  
  return(diff_df)
}

