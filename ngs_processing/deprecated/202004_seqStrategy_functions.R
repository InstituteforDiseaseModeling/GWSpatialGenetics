bc_similarity <- function(df){
  barcode_list <- df2barcodeList(depth_filter(df))
  uniq_pairwise <- data.frame(t(combn(names(barcode_list),2)))
  names(uniq_pairwise) <- c("worm1", "worm2")
  sim_scores <- unlist(lapply(seq(1, length(names(barcode_list))), function(i){
    string_diff(unlist(barcode_list[uniq_pairwise[i, 1]]), unlist(barcode_list[uniq_pairwise[i, 2]]))
  }))
  return(cbind(uniq_pairwise, value=sim_scores))
}

depth_filter <- function(vcf2txt){
  # Subsets SNPs where there are at least 10 reads in each sample 
  dplyr::filter_at(vcf2txt, vars(ends_with(".DP")), all_vars(. > 10)) %>%
    dplyr::select(ends_with(".GT"))
}  

df2barcodeList <- function(filter_df){
  # Turns a dataframe of variants or genotypes into the proper string list format for comparisons.
  # It is presumed samples are the columns mirror VCF format. 
  df_t <- t(filter_df) 
  df_t[df_t == "."] <- "-"
  df_t <- as.data.frame(df_t)
  # concatenate SNP positions into a barcode per sample
  df_t <- tidyr::unite(df_t, barcode, seq(1, ncol(df_t)), sep="")
  barcodes <- as.list(df_t[,1])
  names(barcodes) <- row.names(df_t)
  return(barcodes)
}

gt_count <- function(var_pos){
  var_counts <- cbind.data.frame(
    A=rowSums(var_pos == "A"),
    C=rowSums(var_pos == "C"),
    G=rowSums(var_pos == "G"),
    T=rowSums(var_pos == "T"),
    N=rowSums(var_pos == "."),
    DP_min=rowSums(var_pos < 10))
  var_counts$unique_bases <- rowSums(var_counts[,1:4] != 0)
  var_counts$singleton_bases <- rowSums(var_counts[,1:4] > 1)
  return(cbind.data.frame(var_pos, var_counts))
}

pairwise_plots <- function(sim_df, within_count, outside_count){
  ggplot(sim_df, aes(x=value, group=region, fill=region)) + 
    geom_histogram(aes(y = ..density..), bins=20, alpha=0.2, position="identity") +
    geom_density(alpha=0.5) +
    scale_fill_manual(values = c("pink", "black", "red"),
                      labels = c(paste0("All variants\n(n=", outside_count + within_count, ")"),
                                 paste0("Outside loci\n(n=", outside_count, ")"),
                                 paste0("Within loci\n(n=", within_count, ")")),
                      name = "") + 
    labs(x="Pairwise genetic similarity", y="Comparisons") + 
    theme(legend.position = "top") +
    xlim(c(0,1))
}

sim_scores <- function(pos_df){
  df_t <- tidyr::unite(data.frame(t(pos_df[,c(4:20)])), barcode, seq(1, nrow(pos_df)), sep="")
  var_barcodes <- as.list(df_t[,1])
  names(var_barcodes) <- row.names(df_t)
  uniq_pairwise <- data.frame(t(combn(names(var_barcodes),2)))
  names(uniq_pairwise) <- c("worm1", "worm2")
  sim_scores <- unlist(lapply(seq(1, length(names(var_barcodes))), function(i){
    string_diff(unlist(var_barcodes[uniq_pairwise[i, 1]]), unlist(var_barcodes[uniq_pairwise[i, 2]]))
  }))
  return(cbind(uniq_pairwise, value=sim_scores, bc_length=nchar(var_barcodes[1])))
}

string_diff <- function(a, b){
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  seq.a <- unlist(strsplit(a,split=""))
  seq.b <- unlist(strsplit(b,split=""))
  diff.d <- rbind(seq.a,seq.b)
  only.diff <-diff.d[,diff.d[1,]!=diff.d[2,]]
  pos <- which(diff.d[1,]!=diff.d[2,])
  return(1-(length(pos)/nchar(a)))
}

theme_j <- function () {
  theme_bw(base_size=16) %+replace%
    theme(
      # font sizes and color
      panel.background  = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.background   = element_rect(fill="transparent", colour=NA),
      plot.title        = element_text(size = rel(.85)),
      strip.background  = element_rect(fill="transparent", colour=NA),
      strip.text        = element_text(face="bold", size=rel(.8)),
      axis.title        = element_text(size=rel(0.6)),
      axis.text         = element_text(size=rel(0.5), color="grey30"),
      # legend
      legend.title         = element_text(size=rel(0.8)),
      legend.text          = element_text(size=rel(0.6)),
      legend.background    = element_rect(fill="transparent", colour=NA),
      legend.key           = element_rect(fill="transparent", colour=NA),
      legend.justification = "top"
    )
}

var_group <- function(df){
  tmp_gt <- read.delim(df) %>% dplyr::filter(TYPE == "SNP", !grepl("\\*", ALT)) %>%
    dplyr::mutate(bc_region = ifelse(POS %in% seq(3788, 4534) | 
                                       POS %in% seq(2628, 3345) | POS %in% seq(12562, 14566), "True", "False"))
}