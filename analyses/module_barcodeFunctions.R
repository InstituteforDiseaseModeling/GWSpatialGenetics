#######################################################################
# Functions for barcode analyses 
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: 2025/11
#######################################################################

################################################################################
# Barcode creation and quantifying functions
################################################################################
grpid <- function(x) match(x, unique(x))

gt2barcode <- function(vcf){
  # get exact genotypes 
  tmp_gt <- data.frame(extract.gt(vcf, element = "GT", return.alleles = TRUE))
  print(paste("VCF contains", nrow(tmp_gt), "postions for barcode generation."))
  missing_var = paste0("missing_n_", nrow(tmp_gt))
  
  tmp_gt[is.na(tmp_gt)] <- "X" 
  tmp_gt[tmp_gt == "\\*"] <- "X"
  sequences <- data.frame(sequence = apply(tmp_gt, 2, paste, collapse="")) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::group_by(sequence) %>% 
    dplyr::mutate(!!missing_var := stringr::str_count(sequence, "X"),
                  group = group_indices(),
                  frequency = n(), 
                  analysis_inclusion = ifelse(get(missing_var) >= floor(nrow(tmp_gt) * opt$filt_prop), "Excluded", "Included")) %>%
    dplyr::arrange(desc(frequency), group) %>% ungroup() %>%          
    dplyr::mutate(
      vassar_worm = gsub("^[^_]*_|.batch.*", "", sample),
      amplicon_barcode = group %>% grpid)
  
  return(sequences)
}

################################################################################
# Manual barcode clustering quantification and reassignment
################################################################################
AssignAmpliconGrouping <- function(df, group_variable = "amplicon_barcode"){
  complete_barcodes <- dplyr::select(df, amplicon_barcode_conservative, eval(missing_count_variable)) %>%
    unique() %>%
    filter(eval(missing_count_variable) == 0) %>% .[['amplicon_barcode']]
  
  df_tmp <- df %>%
    dplyr::mutate(confidence = case_when(
      get(missing_count_variable) == 0 ~ "Complete",
      get(missing_count_variable) <= floor(n_variants * 0.01) ~ "High confidence",
      get(missing_count_variable) <= floor(n_variants * 0.02)  ~ "Medium confidence", 
      get(missing_count_variable) >= floor(n_variants * opt$filt_prop) ~ "Excluded",
      TRUE ~ "Lower confidence")) %>% 
    dplyr::group_by(!!! rlang::syms(group_variable)) %>%
    dplyr::mutate(frequency = n(),
      amplicon = case_when(
        confidence == "Excluded" ~ "Excluded",
        (complete_regroup == "Yes" | confidence != "Lower confidence") & frequency >= 10 ~ paste(!!! rlang::syms(group_variable)),
        (complete_regroup == "Yes" | confidence != "Lower confidence") & frequency == 1 ~ "Observed once",
        (complete_regroup == "Yes" | confidence != "Lower confidence") & frequency < 10 ~ "Observed in < 10 specimens",
        confidence == "Lower confidence" & frequency == 1 ~ "Lower confidence - Observed once",
        confidence == "Lower confidence" & frequency > 1 ~ "Lower confidence - Observed > 1"
      ))
  return(df_tmp)
}


BarcodeClusterDifferences <- function(metadata){
  changes <- dplyr::full_join(
    data.frame(table(metadata$amplicon_barcode_conservative)) %>% dplyr::rename(Conservative=Freq),
    data.frame(table(metadata$amplicon_barcode)) %>% dplyr::rename(Clustered=Freq)) %>%
    dplyr::mutate(Var1 = as.numeric(Var1)) %>%
    dplyr::inner_join(.,
      unique(dplyr::select(metadata, amplicon_barcode_conservative, confidence, eval(missing_count_variable))) %>% 
      dplyr::rename(Var1=amplicon_barcode_conservative)) %>%
    dplyr::mutate(Clustered = ifelse(is.na(Clustered), 0, Clustered),
                  difference = Clustered-Conservative,
                  change = case_when(
                    confidence == "Excluded" ~ "Excluded",
                    difference == 0 ~ "No change",
                    difference > 0 ~ "Increased",
                    Clustered == 0 ~ "Regrouped"
                  ),
                  change = factor(change, levels = confidence_change_order),
                  fill_color = ifelse(Conservative == 1 & change != "Excluded", paste0(confidence, " - Observed once"), confidence),
                  fill_color = factor(fill_color, levels = confidence_regroup_order))
  return(changes)
}


################################################################################
# Identification of cross border lineages
# Purpose: Create exploratory groupings for kinship analysis based on shared mitochondrial lineages
################################################################################
Barcode2PotentialKinship <- function(i){
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
}


################################################################################
# Pairwise barcode comparison
################################################################################
# cpp function aided by ChatGPT
cppFunction('
List compare_strings_cpp(CharacterVector strings) {
  int n = strings.size();
  List results;
  
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      std::string s1 = as<std::string>(strings[i]);
      std::string s2 = as<std::string>(strings[j]);
      
      int len = s1.length();
      int matches_with_exclusions = 0;
      int matches_without_exclusions = 0;
      int excluded = 0;
      
      for (int k = 0; k < len; ++k) {
        // Check if position is excluded (one or both characters are "X")
        if (s1[k] == s2[k]){
          matches_without_exclusions++;
        }
        if (s1[k] == s2[k] && !(s1[k] == \'X\' || s2[k] == \'X\')) {
          matches_with_exclusions++;
        } 
        if (s1[k] == \'X\' || s2[k] == \'X\') {
          excluded++;
        }
      }
      
      // Store results for the current pair
      results.push_back(
        List::create(
          Named("pair") = CharacterVector::create(strings[i], strings[j]),
          Named("matches_with_exclusions") = matches_with_exclusions,
          Named("matches_without_exclusions") = matches_without_exclusions,
          Named("excluded") = excluded
        )
      );
    }
  }
  
  return results;
}
')


# R wrapper function
BarcodePairwiseParallel <- function(strings, 
                                    num_cores = detectCores()/2, 
                                    batch_size = 10000) {
  
  cat("Total unique barcodes", length(strings), "\n")
  
  if (!is.character(strings)) {
    stop("Input must be a character vector.")
  }
  
  # Ensure all strings are of the same length
  if (length(unique(nchar(strings))) > 1) {
    stop("All strings must have the same length.")
  }
  
  # Generate all unique pairs of indices
  n <- length(strings)
  all_pairs <- combn(n, 2, simplify = FALSE)
  
  # Split pairs into batches
  batches <- split(all_pairs, ceiling(seq_along(all_pairs) / batch_size))
  # Debug: Verify number of batches
  cat("Total batches created:", length(batches), "\n")
  
  # Function to process a single batch
  process_batch <- function(batch) {
    results <- lapply(batch, function(pair) {
      i <- pair[1]
      j <- pair[2]
      
      result <- compare_strings_cpp(strings[c(i, j)])
      n_variants <- nchar(result[[1]]$pair[1])
      
      # Extract data for the pair
      data.frame(
        sequence.x = result[[1]]$pair[1],
        sequence.y = result[[1]]$pair[2],
        all = 1 - result[[1]]$matches_without_exclusions / n_variants,
        rmNA = 1 - result[[1]]$matches_with_exclusions /(n_variants - result[[1]]$excluded),
        missing = result[[1]]$excluded,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, results)
  }
  
  # Process batches in parallel
  if (.Platform$OS.type == "unix") {
    # Use mclapply for multicore processing on Unix systems
    results <- mclapply(batches, process_batch, mc.cores = num_cores)
  } else {
    # Use parLapply for Windows
    cl <- makeCluster(num_cores)
    results <- parLapply(cl, batches, process_batch)
    stopCluster(cl)
  }
  
  # Combine all batch results into a single data frame
  final_result <- do.call(rbind, results)
  
  return(final_result)
}


################################################################################
# Barcode grouping with missingness
################################################################################
BarcodeAmbiguity <- function(l, pairwise_matrix , match_barcodes = complete_barcodes){
  complete_pairwise_scores <- c(pairwise_matrix[row.names(pairwise_matrix) == l,], pairwise_matrix[,colnames(pairwise_matrix) == l])
  complete_pairwise_scores <- complete_pairwise_scores[!is.na(complete_pairwise_scores) & names(complete_pairwise_scores) %in% match_barcodes]
  potential_match_list <- names(complete_pairwise_scores)[complete_pairwise_scores == 0]
  
  stats = list(amplicon_barcode = l, 
               ambiguous = ifelse(length(potential_match_list) == 0, NA,
                                  ifelse(length(potential_match_list) == 1, "No", "Yes")),
               potential_barcode_matches = case_when(
                 length(potential_match_list) == 0 ~ "None",
                 length(potential_match_list) == 1 ~ potential_match_list[1],
                 TRUE ~ paste(sort(as.numeric(potential_match_list)), collapse=",")
                 )
               )
  return(stats)
}


ManualBarcodeRecode <- function(barcodes, pairs = barcode_pairwise){
  complete_barcodes <- dplyr::filter(barcodes, get(missing_count_variable) == 0) %>% .[['amplicon_barcode']]
  missing_barcodes  <- dplyr::filter(barcodes, get(missing_count_variable) != 0 & 
                                       get(missing_count_variable) < floor(n_variants * opt$filt_prop)) %>% .[['amplicon_barcode']]
  
  # matrix version of relatedness for faster checking of relatedness
  m <- reshape2::acast(pairs, amplicon_barcode.x ~ amplicon_barcode.y, value.var="rmNA")
  
  # Idea - for barcodes with less than 10% of positions missing, check that all pairs in grouping have similarity of 0. 
  # If so, merge as new group to reduce singletons
  refactor_df <- dplyr::bind_rows(lapply(missing_barcodes, function(i) BarcodeAmbiguity(i, m))) 
  
  # check that for any newly grouped barcodes, all barcodes in new group also are identical at available variant positions
  refactor_unambiguous <- dplyr::filter(refactor_df, ambiguous == "No") 
  refactor_counts <- table(refactor_unambiguous$potential_barcode_matches) 
  check_refactor <- names(refactor_counts)[refactor_counts > 1]
  
  exclude_refactor <- unlist(lapply(check_refactor, function(j){
    refactor_tmp <- dplyr::filter(refactor_df, potential_barcode_matches == !!j)
    pairwise_tmp <- m[row.names(m) %in% refactor_tmp$amplicon_barcode,colnames(m) %in% refactor_tmp$amplicon_barcode]
    non_zero_pairs <- sum(pairwise_tmp[pairwise_tmp > 0], na.rm = TRUE)
    if(non_zero_pairs != 0){
      return(j)
    }
  }))
  refactor_unambiguous <- dplyr::filter(refactor_unambiguous, !potential_barcode_matches %in% exclude_refactor)
  
  df <- dplyr::left_join(barcodes, refactor_unambiguous) %>%
    rename(amplicon_barcode_conservative = amplicon_barcode, frequency_conservative = frequency) %>%
    dplyr::mutate(amplicon_barcode = case_when(
      get(missing_count_variable) == 0 | is.na(ambiguous) ~ as.character(amplicon_barcode_conservative),
      ambiguous == "No" ~ as.character(potential_barcode_matches),
      TRUE  ~ "Check"))
  
  rm(m)
  return(df)
}


################################################################################
# Rarefaction estimation
################################################################################
LineageExtrapolation <- function(subset_country, minimum_samples = 30){
  bc_counts <- YearlyBarcodeCounts(subset_country) 
  
  lineage_counts <- dplyr::group_by(bc_counts, year) %>%
    dplyr::summarise(yearly_lineage = n_distinct(amplicon_barcode))
  years_inclusion <- dplyr::filter(lineage_counts, yearly_lineage >= minimum_samples)
  
  yearly_splits <- with(dplyr::filter(bc_counts, year %in% years_inclusion$year), split(barcode_count, year))
  bc_extrapolation <- iNEXT(yearly_splits, q=1, datatype="abundance")
  
  return(bc_extrapolation)
}
