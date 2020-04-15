# guinea worm  - alignment preprocessing functions
#
################################################################################
# creating equilength gene sequences
################################################################################
seq_match <- function(align_dir, align_df, gene){
  seq_df <- data.table::fread(paste(align_dir, align_df, sep="/"))
  names(seq_df) <- c("worm", "start", "sequence")
  seq_df <- dplyr::mutate(seq_df, end = start + nchar(sequence))
  start_min <- min(seq_df$start)
  end_max   <- max(seq_df$end)
  new_sequences <- lapply(setNames(seq_df$worm, seq_df$worm), function(i){
    paste0(paste(rep("-", seq_df[seq_df$worm == i, 2]-start_min), collapse=""), 
           seq_df[seq_df$worm == i, 3], 
           paste(rep("-", end_max - seq_df[seq_df$worm == i, 4]), collapse=""))
  })
  new_sequences <-data.frame(unlist(new_sequences)) %>% tibble::rownames_to_column()
} 


################################################################################
# counting SNPs
################################################################################ 
snp_counter <- function(id_vector, seq_vector, gene, str_pos){
  pos_vector <- gsub("\\-", "N", substring(seq_vector, str_pos, str_pos))
  tmp <- cbind(pos=str_pos, gene=gene, 
               data.frame(table(pos_vector), stringsAsFactors = F))
  names(tmp) <- c("position", "gene", "base", "sample_count")
  if(sum(dplyr::filter(tmp, base != "N")$sample_count == "1") == "1" &
     nrow(dplyr::filter(tmp, base != "N")) == 2){
    alt_base <- dplyr::filter(tmp, base != "N" & sample_count == "1") %>% .[["base"]]
    tmp$singleton_samp = id_vector[pos_vector == alt_base]
  } else{
    tmp$singleton_samp = "None"
  }
  return(tmp)
}

fasta2snpIden <- function(df, gene, start, end, missing_min){
  # no longer reading directly from fasta with potential bug in sam2fasta.py. use input fom seq_match function instead. 
  # fasta    <- Biostrings::readDNAStringSet(file = paste0(path, "/2020_merged_", gene, ".fasta"))
  # fasta_df <- data.frame(worm = names(fasta), sequence = paste(fasta), stringsAsFactors = F)[-1,]
  # count the number of bases at each location
  df_naOmit <- na.omit(df[,c("worm", gene)])
  var_count  <- bind_rows(lapply(seq(start, end), function(j){snp_counter(df_naOmit[,1], df_naOmit[,2], gene, j)}))
  head(var_count)
  # subset columns with at least two variants, excluding missing 
  var_count <- tidyr::spread(var_count, key = base, value=sample_count) %>%
    dplyr::mutate(variant = ifelse((N <= missing_min | is.na(N))  & rowSums(is.na(.[,c("A", "C", "G", "T")])) <= 2, "T", "F"))
}

barcode_create <- function(seq_vector, snp_pos){
  barcode <- character(length(snp_pos))
  for (i in 1:length(snp_pos)){
    barcode[i] <- substring(seq_vector, snp_pos[i], snp_pos[i])
  }
  return(paste(barcode, collapse=""))
}  

sequence2barcode <- function(usr_gene, worm_vector, seq_vector, variant_df){
  filtered_pos <- dplyr::filter(variant_df, gene == usr_gene,  variant == "T" & !singleton_samp %in% sing_outlier)
  singletonInc_bc <- lapply(seq_vector, function(j){barcode_create(j, filtered_pos$position)})
  #singletonExc_bc <- lapply(seq_vector, function(j){barcode_create(j, dplyr::filter(filtered_pos, singleton_samp == "None") %>% .[["position"]])})
  seq_df <- cbind.data.frame(worm=worm_vector, unlist(singletonInc_bc), stringsAsFactors = F)
  names(seq_df) <- c("worm", usr_gene)
  return(seq_df)
}
