################################################################################
# Commonly used functions for Guinea worm genomic analyses
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: May 2023
################################################################################

################################################################################
# General
################################################################################
mutate_cond <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  # Initialize any new variables as new_init
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data %>% filter(condition) %>% mutate(...)
  .data
}


make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
    matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"), 
    interpolate = TRUE
  )
}


################################################################################
# Merge metadata
################################################################################
MergeMeta <- function(df, metadata = metadata_slim){
  df <- df %>%
    dplyr::inner_join(., dplyr::rename(metadata, 'sample.x'='sample'), by="sample.x") %>%
    dplyr::inner_join(., dplyr::rename(metadata, 'sample.y'='sample'), by="sample.y") %>%
    dplyr::mutate(
      country_pair = ifelse(country.x == country.y, country.x, "Different countries"),
      country_pair = factor(country_pair, levels=names(country_colors)),
      host_pair = ifelse(host.x == host.y, host.x, "Different hosts"),
      year_pair = ifelse(year.x == year.y, year.x, "Across years"),
      diff_days = abs(as.numeric(difftime(emergence_date.x, emergence_date.y, unit="days"))),
      time_pair = ifelse(420 >= diff_days & diff_days >= 300, "10-14 Months", "Not in 10-14 months")
    )
  return(df)
} 


################################################################################
# Relatedness functions
################################################################################
HaversineDist <- function(gps_coords){
  d0 <- geodist::geodist(dplyr::select(gps_coords, gps_e, gps_n) %>% 
                           dplyr::rename(longitude = "gps_e", latitude = "gps_n"),
                         measure = "haversine")
  colnames(d0) <- gps_coords$sample
  row.names(d0) <- gps_coords$sample
  return(d0)
}

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
  diff_df$missing = MissingCountsV(na_pos, diff_df$Var1, diff_df$Var2)/ncol(gt_mat)
  diff_df <- dplyr::mutate_if(diff_df, is.numeric, round, digits=3)
  return(diff_df)
} 

LowerMatrix <- function(matrix){
  if(class(matrix) != "matrix"){
    matrix <- as.matrix(matrix)
  }
  matrix[lower.tri(matrix, diag = FALSE)] <- NA
  return(matrix)
}

MissingCounts <- function(pos_list, sample1, sample2){
  pairwise_missing <- unlist(c(pos_list[sample1], pos_list[sample2]))
  return(length(pairwise_missing))
}
MissingCountsV <- Vectorize(MissingCounts, vectorize.args = c("sample1", "sample2"))


################################################################################
# Plotting functions
################################################################################

################################################################################
# colors

# alleles
allele_colors <- setNames(c("#800000FF", "#FFA319FF", "#8A9054FF", "#155F83FF", "#999999", "#000000"),
                          c("A", "C", "G", "T", "Missing", '*'))

# nextstrain barcodes
base_colors <- c("#FF3200FF",  "#D12600", "#DB6A00", "#F9AB0EFF", "#FFED00FF", 
                 "#B2FF2E", "#00AD00", "#019875FF", "#1BB6AFFF", "#32B2DAFF",
                 "#0076BBFF", "#005B94", "#1E2085", "#610052", "#953272", 
                 "#C70E7BFF", "#FC6882FF", "#FF847CFF")
nextstrain_colors <- setNames(c(base_colors, 
                                shades::saturation(base_colors, shades::scalefac(0.50)),
                                shades::saturation(base_colors, shades::scalefac(0.25))),                                
                                seq(1, length(base_colors)*3))
reduced_base <- setNames(c("#D9C6B8", "#C2B0A3", "#836F65", "#52271CFF"), 
                         c("Observed once", "Observed in < 5 samples", "Observed in < 10 samples", "Observed in < 20 samples")) 
extra_colors <- c("#E51E32FF", "#FF782AFF", "#FDA805FF", "#E2CF04FF", "#B1CA05FF", "#98C217FF", "#779815FF", "#029E77FF", "#09989CFF", "#059CCDFF", "#3F64CEFF", "#7E2B8EFF")
extra_colors <- setNames(c(extra_colors, shades::saturation(extra_colors, shades::scalefac(0.50))),
                              seq(length(nextstrain_colors)+1, length(nextstrain_colors)+length(extra_colors)*2))
nextstrain_colors <- c(nextstrain_colors, reduced_base, extra_colors)

# countries
set.seed(15)
# country_colors <- as.vector(paletteer_d("awtools::bpalette"))
country_colors <- c('#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377',
  shades::saturation(c('#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377'), shades::scalefac(0.60)))
countries <- c("Chad", "Ethiopia", "Mali", "South Sudan", "Cameroon", "Angola", "Niger", "Sudan", "Cote d'Ivoire", "Central African Republic", "Burkina Faso", "Ghana")
country_colors <- setNames(c(country_colors, "#666666"),
  c(countries, "Different countries"))


# hosts 
host_colors <- setNames(
  c("#004F7A", "#30B4CC",
    rev(c("#E7E5CCFF", "#C2D6A4FF", "#9CC184FF", "#669D62FF", "#3C7C3DFF", "#1F5B25FF")),
    "#CC3D24", "#F3C558", "#999999", "#666666"),
  c("Human", "Dog", 
    "Cat (domestic)", "Cat (wild unknown spp.)", "Leopard (Panthera pardus)", "Serval", "Genet", "Civet",
    "Baboon (Papio anubis)", "Equine (donkey etc)", 
    "Unknown", "Different hosts"))

  
# c(rev(c("#6DAE90FF", "#CC3D24", "#999999", "#F3C558", "#088158", "#30B4CC", "#004F7A")), "#666666"),
#                       c("Human", "Dog", "Cat, domestic", "Cat, wild unknown spp.", "Baboon (Papio anubis)", "Leopard (Panthera pardus)", "Unknown", "Different hosts"))

################################################################################
# barcode
position_reformat <- function(df, y_var){
  gt_long <- gt[row.names(gt) %in% df$sample,] %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::inner_join(dplyr::select(df, sample, y_var), .) %>%
    tidyr::pivot_longer(cols = starts_with("ENA"), names_to = "position", values_to = "Allele") 
  gt_long$Allele <- gt_long$Allele %>% replace_na('Missing') 
  gt_long$Allele <- factor(gt_long$Allele, levels=c("A", "C", "G", "T", "Missing", '*'))
  return(gt_long)
}


position_plot <- function(df, y_var){
  order <- dplyr::arrange(df, frequency, amplicon_barcode) %>% .[[eval(y_var)]] %>% unique()
  df_long <- position_reformat(df, eval(y_var)) 
  
  df_long %>% 
    ggplot2::ggplot(aes(x=position, 
                        #y=get(y_var),
                        y=factor(get(y_var), levels=order), 
                        fill=Allele)) +
    geom_tile() +
    scale_fill_manual(values = allele_colors) +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank()) +
    labs(x=paste0("Barcode (", length(unique(df_long$position)), " population identified variants)"))
}

################################################################################
# epidemiological data


################################################################################
# relatedness matrices
related_and_missing_plot <- function(df, sample_id="sample", filter_min=0){
  df <- dplyr::mutate(df, id1 = get(paste0(sample_id, ".x")),
                          id2 = get(paste0(sample_id, ".y")))
  if(filter_min != 0){
    df <- dplyr::filter(df, missing <= filter_min)
  }
  
  if(sample_id != "sample"){
    # placeholder for any potential missing information 
    df <- dplyr::mutate(df, 
                       id1= ifelse(is.na(id1), paste("No ID - worm", vassar_worm.x), id1),
                       id2= ifelse(is.na(id2), paste("No ID - worm", vassar_worm.y), id2))
    }
    # keep order of the original sample names, to keep upper/lower matrices in the correct order
    df_match <- rbind(dplyr::select(df, sample.x, id1) %>% rename(sample = sample.x, id = id1),
                      dplyr::select(df, sample.y, id2) %>% rename(sample = sample.y, id = id2)) %>%
      unique() %>%
      dplyr::arrange(sample) %>%
      dplyr::mutate(id = as.factor(id))

  tmp_pair_related <- dplyr::select(df, id1, id2, rmNA) %>%
    dplyr::rename("sample1" = id1, "sample2" = id2) %>%
    dplyr::mutate(value = 1 - rmNA, 
                  sample1 = factor(sample1, levels=df_match$id),
                  sample2 = factor(sample2, levels=df_match$id))
  # reverse order to fill in bottom matrix
  tmp_pair_missing <- dplyr::select(df, id1, id2, missing) %>%
    dplyr::rename("sample1" = id2, "sample2" = id1) %>%
    dplyr::mutate(value = 1 - missing,
                  sample1 = factor(sample1, levels=(df_match$id)),
                  sample2 = factor(sample2, levels=(df_match$id)))
  
  p_matrix <- tmp_pair_related %>%
    ggplot(aes(x=sample1, y=sample2)) +
    geom_tile(aes(fill=value)) +
    geom_tile(data=dplyr::filter(tmp_pair_related, rmNA == 0),
              aes(x=sample1, y=sample2),
              fill="transparent", colour="black", size=1) +
    geom_text(aes(label = round(value, 3)), size=2) +
    scale_fill_viridis_c(limits = c(0.75, 1), option = 'D', 
                         name="Mitochondrial\ngenetic similarity\n(Excluding missing\npositions)") +
    
    # add the pairwise completeness
    ggnewscale::new_scale_fill() + 
    geom_tile(data = tmp_pair_missing, aes(fill = value)) +
    geom_text(data = tmp_pair_missing, aes(label = round(value, 3)), size=2) +
    scale_fill_viridis_c(limits = c(0.5, 1), option="magma", 
                         name="Pairwise\ncompleteness\t\t")
  
  clean_matrix <- p_matrix +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 0)) + 
    scale_x_discrete(expand = c(0, 0),position = 'top') +
    coord_fixed()
  
  return(clean_matrix)
}

# relatedness distributions
similarity_distributions <- function(df, similarity = "rmNA"){
  x_axis <- ifelse(similarity == "All", 
                   "mtDNA genetic similarity\n(Including missing positions for pairs of worms)", 
                   "mtDNA genetic similarity\n(Excluding missing positions for pairs of worms)")
  
  p_dist <- df %>%
    ggplot(aes(x=1-get(similarity), y= ..scaled.., fill=country_pair)) + 
    geom_density(alpha=0.5) +
    scale_fill_manual(values=country_colors) +
    labs(x=x_axis, y="Scaled probability density", 
         fill="Country of pairs of worms")
  
  g <- make_gradient(deg = 45, n = 200, cols = RColorBrewer::brewer.pal(9, "Greys")[1:5])
  c_dist <- df %>% 
    ggplot(aes(x=1-get(similarity))) + 
    annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    stat_ecdf(aes(color = country_pair)) +
    scale_color_manual(values=country_colors) +
    labs(x=x_axis, y="Cumulative density",
         color="Country of pairs of worms") +
    annotate(geom="text", label="Pairs are\nless related", x=0.75, y=0.95, color="black") +
    annotate(geom="text", label="Pairs are\nmore related", x=0.95, y=0.10, color="black")
  
  density_p <- ggpubr::ggarrange(p_dist, c_dist, ncol=2, common.legend = T)
  return(density_p)
} 

# spatial relatedness
dist_by_sim_plot <- function(df, diff_all = genetic_diff){
  df <- df  %>%
    dplyr::mutate(year = as.character(year),
                  has_gps = ifelse(!is.na(gps_n) & !is.na(gps_e), "Provided", "Missing"))
  sample_counts <- df %>% ggplot(aes(x=year, y=..count.., fill=host)) +
    geom_bar(color="black", size=0.5) +
    scale_fill_manual(name="", values=host_colors) +
    labs(title = "Species of host", x="Year", y="Sequenced specimens") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  coord_counts <- df %>% ggplot(aes(x=year, y=..count.., fill=has_gps)) +
    geom_bar(color="black", size=0.5) +
    scale_fill_manual(name="", values=c("ivory3", "seagreen3")) +
    labs(title = "Specimens with GPS coordinates", x="Year", y="Sequenced specimens") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  bins <- dplyr::filter(diff_all, sample.x %in% df$sample & sample.y %in% df$sample) %>%
    MergeMeta(.) %>%
    left_join(., dplyr::rename(pairwise_dist, sample.x=Var1, sample.y=Var2)) %>% 
    dplyr::mutate(year_diff=abs(as.numeric(year.x) - as.numeric(year.y)),
                  year_diff=ifelse(year_diff == 0, "Same year",
                                   ifelse(year_diff == 1, "Within one year", "> One year"))) %>%
    ggplot(aes(x=meters/1000,y=1-rmNA)) + 
    geom_point(aes(color=as.character(year_diff)), alpha=0.75) +
    scale_colour_viridis_d(option = "inferno") +
    #stat_binhex() +
    # scale_fill_viridis_c(trans = "log", option = "inferno", na.value = "grey50") +
    labs(x="Distance between pairs of worm samples\n(kilometers)",
         y="Pairwise mtDNA genetic similarity\n(Excluding missing positions in barcode)",
         fill="Number of pairwise\ncomparisons",
         color="Years difference\nbetween pairs",
         title="Association between pairwise distance\nand genetic similarity for pairs of worms") + 
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 30, hjust = 1))
  p <- ggpubr::ggarrange(
    ggpubr::ggarrange(sample_counts, coord_counts, nrow=2, align = "hv"), 
    bins,
    nrow=1
  )
  return(p)
}  

################################################################################
# multipanel plots
ClusterPlots <- function(df, cluster_id=NA, output_name=NA, filter_min=0, 
                         pairwise_df = genetic_diff,
                         #pairwise_df = diff_all, 
                         sample_name = "sample"
                         #sample_name="wormnum_dpdxid"
                         ){
  if(is.na(output_name)){
    if(!is.na(cluster_id)){
      output_name <-  paste0("cluster_", gsub(" |cluster", "", cluster_id))
    } else{
      output_name <- "unnamedSimilarity"
    }
  }
  
  print(paste("Checking", output_name))
  # subset sampled in cluster
  if(!is.na(cluster_id)){
    df_cluster <- dplyr::filter(df, cluster == !!cluster_id)
  } else {
    df_cluster <- df
  }
  
  if(nrow(df_cluster) > 2 & sum(!is.na(df_cluster$amplicon)) > 2){
    # general metrics
    host <- df_cluster %>% ggplot(aes(x=as.character(year), y=..count.., fill=host))+
      geom_bar(colour="black") +
      scale_fill_manual(values = host_colors, name="Host") +
      labs(x="Year", y="Specimens") + 
      theme(legend.position="top") +
      ggtitle("")
    barcodes <- df_cluster %>% ggplot(aes(x=as.character(year), y=..count.., fill=amplicon))+
      geom_bar(colour="black") +
      scale_fill_manual(values = nextstrain_colors,
                        name="Nextstrain\ngroup") +
      labs(x="Year", y="Specimens") + 
      theme(legend.position="top") +
      guides(fill=guide_legend(nrow=2, byrow=TRUE))
    
    general_joint <- ggpubr::ggarrange(host, barcodes, nrow = 2, widths = c(1,1), align="hv")
    
    # relatedness and missing-ness plots
    tmp_pair <- dplyr::filter(pairwise_df, sample.x %in% df_cluster$sample & sample.y %in% df_cluster$sample) %>%
      mutate(outline = ifelse(rmNA == 0, T, NA)) %>%
      MergeMeta(.)
    related_joint <- related_and_missing_plot(tmp_pair, eval(sample_name), filter_min=filter_min)
    
    # variant position plots
    # position <- position_plot(df_cluster, eval(sample_name))
    
    # all combined
    plot_name <- paste0("Cluster ID: ", output_name, "\n", 
                        "Narrative samples: ", nrow(df_cluster), "; Sequenced specimens: ", sum(!is.na(df_cluster$amplicon)))
    plot <-  ggpubr::ggarrange(
      ggpubr::ggarrange(general_joint, related_joint, ncol=2, widths = c(1,2))) #,
      # position, nrow=2, heights = c(2.75,1))
    plot <- ggpubr::annotate_figure(plot, top = plot_name)
    
    # save plot
    plot_height <- ifelse(nrow(df_cluster) <= 15, 8, 21) 
    plot_width <- ifelse(nrow(df_cluster) <= 15, 10, 26)
    SavePlots(plot, output_dir, paste0("relatednessMatrix_", output_name, ".png"), height=plot_height, width=plot_width)
    
    potential_outliers <- "Exit"
  } else{
    potential_outliers <- "Skipped"
  }
  return(potential_outliers)
} 


################################################################################
# saving
SavePlots <- function(ggplot_obj, plot_path, name, width=8, height=4){
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name), 
         plot = ggplot_obj,
         path = plot_path,
         width = width, height = height, units = c("in"), dpi = 100)
}