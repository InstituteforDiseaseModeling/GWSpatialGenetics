################################################################################
# Commonly used functions for Guinea worm genomic analyses
# Plotting functions
# Author: Jessica Ribado - Institute for Disease Modeling 
# Date: May 2023, updated February 2025
################################################################################

################################################################################
# Misc
################################################################################
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
# Color schemes
################################################################################
# alleles
allele_colors <- setNames(c("#800000FF", "#FFA319FF", "#8A9054FF", "#155F83FF", "#999999", "#000000"),
                          c("A", "C", "G", "T", "Missing", '*'))

# nextstrain barcodes
base_colors <- c("#FF3200FF",  "#D12600", "#DB6A00", "#F9AB0EFF", "#FFED00FF", 
                 "#B2FF2E", "#00AD00", "#019875FF", "#1BB6AFFF", "#32B2DAFF",
                 "#0076BBFF", "#005B94", "#1E2085", "#610052", "#953272", 
                 "#C70E7BFF", "#FC6882FF", "#FF847CFF")
base_colors2 <- c(base_colors, 
  shades::saturation(base_colors, shades::scalefac(0.50)),
  shades::addmix(base_colors, "blue", 0.5),
  #shades::saturation(shades::addmix(base_colors, "blue", 0.5), shades::scalefac(0.50)),
  shades::addmix(base_colors, "yellow", 0.2),
  shades::saturation(shades::addmix(base_colors, "yellow", 0.2), shades::scalefac(0.50)),
  "#0054FFFF", "#3299FFFF", "#65CCFFFF", "#99EDFFFF", "#CCFFFFFF", "#FFFFCCFF", "#FFEE99FF", "#FFCC65FF", "#FF9932FF", "#FF5500FF",
  "#E76254FF", "#EF8A47FF", "#F7AA58FF", "#FFD06FFF", "#FFE6B7FF", "#AADCE0FF", "#72BCD5FF", "#528FADFF", "#376795FF", "#1E466EFF")
all_colors <- setNames(base_colors2, seq(1, length(base_colors2)))


# reduced_base <- setNames(c("#D9C6B8", "#C2B0A3", "#836F65", "#52271CFF"), 
#                          c("Observed once", "Observed in < 5 samples", "Observed in < 10 samples", "Observed in < 20 samples")) 
reduced_base <- setNames(
  c("#836F65", "#D9C6B8", "#B7957CFF", "#8C6751FF", "#593527FF", "#888888", "#444444"), 
  c("Observed in < 20 specimens", "Observed in < 10 specimens", "Observed once", 
    "Lower confidence - Observed > 1", "Lower confidence - Observed once", "Excluded", "Not sequenced")) 
nextstrain_colors <- c(all_colors, reduced_base)

# countries
set.seed(15)
# country_colors <- as.vector(paletteer_d("awtools::bpalette"))
country_colors <- c('#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377',
                    shades::saturation(c('#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377'), shades::scalefac(0.60)))
countries <- c("Chad", "Ethiopia", "Mali", "South Sudan", "Cameroon", "Angola", "Niger", "Sudan", "Cote d'Ivoire", "Central African Republic", "Burkina Faso", "Ghana")
country_colors <- setNames(c(country_colors, "#666666", "#666666"),
                           c(countries, "Different countries", "Mixed"))


# hosts 
host_colors <- setNames(c("#004F7A", "#30B4CC", "#BFEFFF", "#F3C558", "#FFA500",
                          "#C2D6A4FF", "#9CC184FF", "#669D62FF", "#3C7C3DFF", "#1F5B25FF", "#1E3D14FF", "#192813FF",
                          "#C0431FFF", "#A0522D",
                          "#FAEBD7",
                          "#999999", "#666666", "#666666"),
                        c("Human", "Baboon (Papio Anubis)", "Primate (Non-Baboon Non-Human Species)", "Dog", "Wild Canid (Jackal Etc.)",
                          "Cat (Domestic)", "Cat (Wild Unknown Spp.)", "Civet", "Genet", "Leopard (Panthera Pardus)", "Serval", "African Wild Cat",
                          "Equine (Donkey Etc)", "Antelopinae (Gazelles Antelopes Etc)",
                          "Lab Ferret",
                          "Unknown", "Different hosts", "Mixed"))

       

################################################################################
# Barcode quality plots
################################################################################
AlleleFrequencyPlot <- function(input_gt, exclusion_samples){
  gt <- extract.gt(input_gt, element = "GT", return.alleles = TRUE)
  print(paste("Size of VCF pre filtering:", dim(gt)))
  gt[gt == "*"] <- NA
  gt <- gt[,!colnames(gt) %in% exclusion_samples]
  print(paste("Size of VCF post filtering:", dim(gt)))
  
  # get allele counts
  gt_count <- reshape2::dcast(reshape2::melt(gt), Var1 ~ value) %>%
    mutate(across(-1)/rowSums(across(-1))) %>%
    rowwise() %>%
    mutate(Major=max(c_across(A:`T`)),
           Minor=1-Major-`NA`) %>%
    dplyr::rename("Missing" = `NA`) %>%
    dplyr::arrange(desc(Minor), Missing)
  
  # plot
  af_p <- gt_count %>% 
    dplyr::select(Var1, starts_with("M")) %>%
    tidyr::pivot_longer(cols = -Var1, names_to = 'Allele', values_to = 'Proportion') %>%
    dplyr::mutate(Var1 = factor(Var1, levels = gt_count$Var1)) %>%
    ggplot(aes(x=Var1, y=Proportion, fill=Allele)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("grey80", "grey60", "black")) +
    labs(title = paste0("Allele frequencies for ", strsplit(batch_name, "_")[[1]][2], " sequencing batch"),
         x = paste0("Variant positions (n = ", nrow(gt), ")"),
         y = paste0("Proportion of samples (n = ", ncol(gt), ")")) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  return(af_p)
}


MissingFrequencyPlot <- function(barcodes, threshold = opt$filt_prop){
  if(!"missing_prop" %in% names(barcodes)){
    missing_count_variable = names(barcodes)[grepl("missing_n_", names(barcodes))]
    n_variants <- readr::parse_number(missing_count_variable)
    barcodes <- dplyr::mutate(barcodes, missing_prop = get(missing_count_variable)/n_variants)
  }
  
  missing_p <- barcodes %>%
    ggplot(aes(x=missing_prop, y=Conservative)) +
    geom_jitter(alpha=0.5) +
    scale_y_log10() +
    geom_vline(xintercept = threshold, linetype="dashed", color = "red") +
    labs(title = paste0(gsub("batch_", "Specimen library through ", batch_name), " sequencing batch"), 
         x = paste0("Proportion of missing variant calls in the barcode\n(", n_variants, " barcode variants)"),
         y = paste0("Frequency of barcode in sequenced library\n(", sum(barcodes$Conservative, na.rm=T), " specimens)")) +
    annotate(geom = "text", color = "red",
             label = paste0("Sample missingness\nexclusion threshold\n(", sum(barcodes$missing_prop > threshold, na.rm=T), " specimens excluded)"),
             x = threshold, y = 100, hjust = -0.1)
  
  return(missing_p)
}


confidence_bin_subtitle <- paste0("High (< 1%), Medium (< 2%), Lower (<", opt$filt_prop*100, "%), Excluded (>", opt$filt_prop*100, "%)")
confidence_order = c("Complete", "High confidence", "Medium confidence", "Lower confidence", "Excluded")
confidence_change_order  = c("No change", "Increased", "Regrouped", "Excluded")
confidence_regroup_order = c("Complete", "Complete - Observed once", 
                             "High confidence", "High confidence - Observed once",
                             "Medium confidence", "Medium confidence - Observed once", 
                             "Lower confidence", "Lower confidence - Observed once", "Excluded")
confidence_group_colors <- setNames(c(#"#004949FF", "#009292FF", 
                                      "#006DDBFF", "#6DB6FFFF",
                                      "#490092FF", "#B66DFFFF", 
                                      "#FF6DB6FF", "#FFB6DBFF",
                                      "#593527FF", "#8C6751FF",
                                      "#888888"), confidence_regroup_order)

BarcodeConfidencePlot <- function(barcodes){
  barcode_conf_p <- barcodes %>%
    tidyr::pivot_longer(cols = c("Conservative", "Clustered"), names_to = "group", values_to = "Frequency") %>%
    dplyr::mutate(group = factor( group, levels = c("Conservative", "Clustered"))) %>%
    dplyr::filter(Frequency > 0) %>%
    #dplyr::mutate(confidence = factor(confidence, levels = confidence_order)) %>%
    ggplot(aes(x=group, y=Frequency, fill=fill_color)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = confidence_group_colors) +
    labs(x="",
         y=paste0("Barcode frequency (n=", sum(barcodes$Conservative), " specimens)"),
         fill="Barcode\nconfidence\nbins")
  return(barcode_conf_p)
}


BarcodeRegroupChangesPlot <- function(barcodes){
  regroup_changes_p <- barcodes %>% 
    ggplot(aes(x=change, y=after_stat(count), fill = fill_color)) +
    geom_bar(color="black") +
    scale_fill_manual(values=confidence_group_colors) +
    labs(x="Clustering changes", 
         y=paste0("Original barcodes (n=", length(unique(bc_changes$Var1)), ")"),
         fill="Pre-regrouping\nconfidence category")
  return(regroup_changes_p)
}


UniqueBarcodeChangesPlot <- function(barcodes){
  barcode_count_p <- data.frame(Group = factor(c("Conservative", "Clustered"), levels = c("Conservative", "Clustered")),
                                   Counts = c(sum(barcodes$Conservative > 0), sum(barcodes$Clustered > 0))) %>%
    ggplot(aes(x=Group, y=Counts)) +
    geom_bar(stat = "identity") +
    geom_text(stat='identity', aes(label=Counts), vjust=-1) +
    labs(x="", y="Unique barcode groups")
  return(barcode_count_p)
}


BarcodeChangeMulti <- function(barcodes, n_barcode_variants = n_variants){
  confidence_bin_title <- paste0("Barcode confidence based on completeness (n=", n_barcode_variants, " variant positions)")
  plotA <- UniqueBarcodeChangesPlot(barcodes) + ggtitle("Total number of barcodes\ndenoting maternal lineages")
  plotB <- ggpubr::ggarrange(BarcodeConfidencePlot(barcodes) + 
                               guides(fill = guide_legend(nrow = 4)) + 
                               ggtitle("Quality of barcodes for all\nsequenced specimens"),
                             BarcodeRegroupChangesPlot(barcodes) + 
                               ggtitle("Change of barcodes\nfrom conservative categories"),
                             ncol = 2, common.legend = T, widths = c(1.5,2))
  plot <-  ggpubr::ggarrange(plotA, plotB, nrow=1, widths = c(1,2))
  plot <- ggpubr::annotate_figure(plot, top = paste0(confidence_bin_title, "\n", confidence_bin_subtitle))
  return(plot)
}


################################################################################
# Barcode frequency plots
################################################################################
BarcodeByCountryPlot <- function(metadata){
  specimen_counts <- table(metadata$country)
  high_specimens_countries <- names(specimen_counts)[specimen_counts > 10]
  
  barcode_counts <- metadata %>%
    dplyr::filter(!is.na(amplicon) & analysis_inclusion == "Included" & country %in% high_specimens_countries) %>%
    ggplot(aes(x=amplicon, y=after_stat(count), fill=country)) +
    scale_fill_manual(values=country_colors) +
    geom_bar(color="black", linewidth=0.25) +
    labs(x="Barcode group", y="Number of specimens", fill="Country") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  return(barcode_counts)
}


PlotTitle <- function(df){
  countries <- unique(df$country)
  if(length(countries == 1)){
    admins <- unique(df$admina)
    if(length(admins) == 1){
      ggplot_title = paste(gsub(" ", "", countries[1]), admins[1], sep='-')
    } else {
      ggplot_title = countries[1]
    }
  }
  return(ggplot_title)
}

DynamicPlotDimensions <- function(column, num_entries, num_facets){
  if(column == "year"){
    base_height <- 5  # Default height per facet
    base_width  <- 7   # Default width
  } else{
    base_height <- 3.5  # Default height per facet
    base_width  <- 9 
  }

  if(num_entries > 1000 & num_facets == 1){
    base_height <- base_height + 1.5
    base_width  <- base_width + 2.5
  }
  
  if(num_facets > 1){
    base_height <- base_height*ceiling(sqrt(num_facets))/2.5
    base_width <- base_width*ceiling(sqrt(num_facets))/2.5
  }
  return(list(height=base_height, width=base_width))
}


BarcodeCountPlots <- function(df, column = "year", 
                              ggplot_title = NULL, 
                              facet_variable = NULL, 
                              input_width=NULL, input_height=NULL){
    
  base_p <- function(p_df = df, column, facet = facet_variable){
    title_base = "Maternal lineages by barcode"
    batch_reformat = lubridate::mdy(gsub("batch_", "", batch_name))
    
    p <- p_df %>%
      dplyr::filter(!is.na(!!column)) %>%
      ggplot(aes(x=!!column, y=after_stat(count), fill=amplicon)) +
      geom_bar(color="black", linewidth=0.25) +
      scale_fill_manual(values=nextstrain_colors) +
      labs(title = ifelse(is.null(ggplot_title), title_base, paste(title_base, "in", ggplot_title)),
           subtitle = paste("Barcode groups determined by specimens included through", batch_reformat, "sequencing batch"),
           y="Specimens", fill="Maternal lineage\nbarcodes") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    # Conditionally add facet_grid if facet_variable is not NULL
    if (!is.null(facet)) {
      unique_facets <- length(unique(p_df[[facet]]))
      p <- p + facet_wrap(vars(!!sym(facet)),
                          ncol = ifelse(column == "year", ceiling(sqrt(unique_facets)), 1),
                          scales = "free_y")
    }
    return(p)
  }
  
  ggplot_title = ifelse(is.null(ggplot_title), PlotTitle(df), ggplot_title)
  
  if(column != "year"){
    print("Checking the emergence date to plot maternal lineage plots by month.")
    column = "emergence_date_month"
    # Convert all dates to first of the month to group cases by month 
    df[[column]] = date(format(df$emergence_date, "%Y-%m-01"))
    p <- base_p(column = sym(column)) +
      scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
      labs(x="Emergence month and year")
  } else {
    time_seq <- seq(min(df[[column]], na.rm = T), max(df[[column]], na.rm = T))
    p <- base_p(column = sym(column)) +
      scale_x_continuous(breaks = c(time_seq)) +
      labs(x="Emergence year")
  }
  
  num_facets <- if (!is.null(facet_variable)) dplyr::n_distinct(df[[facet_variable]]) else 1
  if(is.null(input_width) | is.null(input_height)){
    output_size <- DynamicPlotDimensions(sym(column), nrow(df), num_facets)
  } else{
    output_size <- list()
    output_size$width = input_width
    output_size$height = input_height
  }
  output_name <- ifelse(is.null(ggplot_title), paste("Unnamed", format(Sys.time(), "%H:%M:%S"), sep="-"), ggplot_title)
  output_name = paste(gsub(" ", "", output_name), gsub("emergence_date_", "", column), sep="_")
  SavePlots(p, output_dir, paste(output_name, "barcodeCounts.png", sep="_"),
           height=output_size$height, width=output_size$width)
  rm(p)
} 
  

BarcodeProportionPlots <- function(p_df){

  count_p <- p_df %>%
    ggplot(aes(x=column, y=..count.., fill=host)) +
    geom_bar(color="black", size=0.25) +
    scale_fill_manual(values=host_colors) +
    scale_x_continuous(breaks = c(seq(min(p_df$year), max(p_df$year)))) +
    labs(title = paste("Specimens from", i),
         x=column, y="Specimens", fill="Host")
   

 prop_p <- p_df %>%
   ggplot(aes(x=year, y=..count.., fill=amplicon)) +
   geom_bar(position = "fill", color="black", size=0.25) +
   scale_fill_manual(values=nextstrain_colors) +
   scale_x_continuous(breaks = c(seq(min(p_df$year), max(p_df$year)))) +
   labs(title = paste("Barcodes in", i),
        subtitle = "(Groups determined by frequency in an all country library)",
        x="Year", y="Proportion of barcodes", fill="Barcode") +
   theme(axis.text.x=element_text(angle=45, hjust=1))
 
 prop_comb <-  ggpubr::ggarrange(
   count_p + theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank()),
   prop_p, nrow=2, heights = c(.5,1), align = "v")
 SavePlots(prop_comb, output_dir, paste0(gsub(" ", "", i), "_year_barcodeProportion.png"), 
           height=8, width=6)
}


################################################################################
# Related matrices
################################################################################
related_and_missing_plot <- function(df, sample_id="sample", 
                                     filter_min = 0,
                                     n_barcode_variants = n_variants){
  
  # set the order by name of the user specified variable for plotting
  df <- df %>%
    rowwise() %>%
    dplyr::mutate(
      id1 = ifelse(get(paste0(sample_id, ".x")) < get(paste0(sample_id, ".y")), get(paste0(sample_id, ".x")), get(paste0(sample_id, ".y"))),
      id2 = ifelse(get(paste0(sample_id, ".x")) >= get(paste0(sample_id, ".y")), get(paste0(sample_id, ".x")), get(paste0(sample_id, ".y"))),
      ) %>%
    dplyr::arrange(id1, id2)
  if(max(df$missing, na.rm = T) > 1){
    df <- dplyr::mutate(df, missing = missing/n_barcode_variants)
  }
  
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
  tmp_pair_related <- dplyr::select(df, id1, id2, rmNA) %>%
    dplyr::rename("sample1" = id1, "sample2" = id2) %>%
    dplyr::mutate(value = 1 - rmNA)
  # reverse order to fill in bottom matrix
  tmp_pair_missing <- dplyr::select(df, id1, id2, missing) %>%
    dplyr::rename("sample1" = id2, "sample2" = id1) %>%
    dplyr::mutate(value = ifelse(sample1 == sample2, NA, 1 - missing))
  
  p_matrix <- tmp_pair_related %>%
    ggplot(aes(x=sample1, y=sample2)) +
    geom_tile(aes(fill=value)) +
    geom_tile(aes(x=sample1, y=sample2), 
              fill="transparent", color="#666666", size=0.1) +
    geom_tile(data=dplyr::filter(tmp_pair_related, rmNA == 0),
              aes(x=sample1, y=sample2),
              fill="transparent", colour="black", size=1) +
    #geom_text(aes(label = round(value, 3)), size=2) +
    geom_text(data = dplyr::filter(tmp_pair_related, value == 1),
              aes(label=value), size=2) +
    scale_fill_viridis_c(limits = c(0.75, 1), option = 'D', 
                         name="Mitochondrial\ngenetic similarity\n(Excluding missing\npositions)") +
    
    # add the pairwise completeness
    ggnewscale::new_scale_fill() + 
    geom_tile(data = tmp_pair_missing, aes(fill = value)) +
    geom_tile(aes(y=sample1, x=sample2), 
              fill="transparent", color="#666666", size=0.1) +
    # geom_text(data = tmp_pair_missing, aes(label = round(value, 3)), size=2) +
    geom_text(data = dplyr::filter(tmp_pair_missing, value == 1),
              aes(label = value), size=2) +  
    scale_fill_viridis_c(limits = c(0.5, 1), option="magma", 
                         name=paste0("Pairwise\ncompleteness\t\t\n(n=", n_barcode_variants, " variants)")) 
  
  clean_matrix <- p_matrix +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0)) + 
    coord_fixed()
  
  return(clean_matrix)
}

# relatedness distributions
similarity_distributions <- function(df, similarity = "rmNA"){
  x_axis <- ifelse(similarity == "All", 
                   "mtDNA genetic similarity\n(Including missing positions for pairs of worms)", 
                   "mtDNA genetic similarity\n(Excluding missing positions for pairs of worms)")
  
  p_dist <- df %>%
    ggplot(aes(x=1-get(similarity), weight=pairwise_count/sum(pairwise_count), fill=country_pair)) + 
    geom_density(alpha=0.5) +
    scale_fill_manual(values=country_colors) +
    labs(x=x_axis, y="Scaled probability density", 
         fill="Country of pairs of worms")

  # Compute cumulative probability manually
  df <- df %>%
    group_by(country_pair) %>%
    arrange(desc(get(similarity))) %>%  # Ensure data is sorted within groups
    mutate(
      cum_count = cumsum(pairwise_count),  # Cumulative sum per group
      cdf = cum_count / sum(pairwise_count)  # Normalize to get cumulative probability
    ) %>% ungroup()
  
  g <- make_gradient(deg = 45, n = 200, cols = RColorBrewer::brewer.pal(9, "Greys")[1:5])
  c_dist <- df %>%
    ggplot(aes(x=1-get(similarity), y=cdf, color=country_pair)) +
    annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    # stat_ecdf(aes(color = country_pair)) +
    geom_step(linewidth = 1) +  
    scale_color_manual(values=country_colors) +
    labs(x=x_axis, y="Cumulative density",
         color="Country of pairs of worms") +
    annotate(geom="text", label="Pairs are\nless related", x=0.85, y=0.95, color="black") +
    annotate(geom="text", label="Pairs are\nmore related", x=0.98, y=0.10, color="black")

  density_p <- ggpubr::ggarrange(p_dist, c_dist, ncol=2, common.legend = T)
  return(density_p)
} 


# spatial relatedness
dist_by_sim_plot <- function(df, pairwise_df){
  df <- df %>%
    dplyr::mutate(year = as.character(year),
                  has_gps = ifelse(!is.na(gps_n) & !is.na(gps_e), "Provided", "Missing"))
  sample_counts <- df %>% ggplot(aes(x=year, y=after_stat(count), fill=host)) +
    geom_bar(color="black", linewidth=0.5) +
    scale_fill_manual(name="", values=host_colors) +
    labs(title = "Species of host", x="Year", y="Sequenced specimens") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  coord_counts <- df %>% ggplot(aes(x=year, y=after_stat(count), fill=has_gps)) +
    geom_bar(color="black", linewidth=0.5) +
    scale_fill_manual(name="", values=setNames(c("ivory3", "seagreen3"), c("Missing", "Provided"))) +
    labs(title = "Specimens with GPS coordinates", x="Year", y="Sequenced specimens") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  bins <- dplyr::filter(pairwise_df, sample.x %in% df$sample & sample.y %in% df$sample) %>%
    dplyr::mutate(year_diff = ifelse(!grepl("_", year_pair), "Same year", "Different years")) %>%
    ggplot(aes(x=meters/1000, y=1-rmNA)) + 
    geom_point(aes(color=as.character(year_diff)), alpha=0.75) +
    scale_colour_manual(values= c("#FFA400FF", "#862633FF", "black")) +
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
    bins, nrow=1)
  
  return(p)
}  


similarity_bin_order <- sort(1 - c(seq(0, 0.1, by=0.01), 0.15, 0.20, 0.25))
PairwiseProportionHeatMap <- function(df){
  # summarize counts
  summary <- df %>%
    dplyr::mutate(similarity_bin = cut(1-rmNA, 
                                       breaks = similarity_bin_order, 
                                       labels = head(similarity_bin_order, -1)), 
                  similarity_bin = ifelse(rmNA == 0, "1", paste(similarity_bin))) %>%
    dplyr::group_by(year_pair, similarity_bin) %>%
    dplyr::summarize(count = sum(pairwise_count)) %>%
    dplyr::group_by(year_pair) %>%
    dplyr::mutate(proportion=count/sum(count)*100) 
  
  p <- summary %>%
    ggplot(aes(x=year_pair, y=similarity_bin, fill=proportion))+
    geom_tile() +
    scale_fill_gradientn(colours=c("#1D3141FF", "#096168FF", "#209478FF", "#75C56EFF", "#E2EE5EFF"),
                         na.value = "transparent",
                         limits=c(0,100),
                         guide = guide_colourbar(
                           title = "Proportion of pairs\nper year",
                           title.position="top", title.hjust = 0.5)) +
    geom_text(aes(label=count)) +
    labs(x="Year", 
         y="Pairwise mtDNA genetic similarity bins\n(Excluding missing positions in barcode)") +
    new_scale_fill() + 
    geom_tile(
      data = data.frame(year_pair = "",
                        similarity_bin = factor(similarity_bin_order, levels = similarity_bin_order),
                        proportion = similarity_bin_order), 
      aes(x=year_pair, y=similarity_bin, fill=proportion)) +
    scale_fill_gradientn(colours=c("#1E313EFF", "#4E475FFF", "#8B5975FF", "#C86C7CFF", "#FA8975FF"),
                         breaks=c(1,0.75), 
                         labels=c("Potentially\t\nrelated\t", "Not\nrelated"),
                         limits=c(0.75,1),
                         guide = guide_colourbar(
                           title = "Pairwise\nrelatedness",
                           title.position="top", title.hjust = 0.5)) +
    theme(legend.position = "top")
  return(p)
}


AddPairwiseData <- function(df, dt, variable){
  lookup_vector <- setNames(dplyr::pull(df, {{variable}}), dplyr::pull(df, sample))
  
  dt <- dt[, (paste(variable, "pair", sep="_")) := data.table::fcase(
    lookup_vector[match(sample.x, names(lookup_vector))] == lookup_vector[match(sample.y, names(lookup_vector))], as.character(lookup_vector[match(sample.x, names(lookup_vector))]),
    is.na(lookup_vector[match(sample.x, names(lookup_vector))]) | is.na(lookup_vector[match(sample.y, names(lookup_vector))]), "Missing data",
    default = "Mixed")]
  return(dt)
}


VariableDistributionPlot <- function(df, dt, variable,
                                     variable_title = NULL,
                                     color_vector = NULL){
  if(variable != "year"){
    variable_order <- names(sort(table(dplyr::pull(df, {{variable}})), decreasing = T))
    df[[variable]] <- factor(df[[variable]], levels=variable_order)
  }
  count_p <- ggplot(df, aes(y=!!sym(variable), x=after_stat(count))) +
    geom_bar() +
    scale_y_discrete(limits=rev) +
    labs(y="", x="Specimen count")
  
  pair_variable <- paste(variable, "pair", sep="_")
  if(!pair_variable %in% names(dt)){
    dt <- AddPairwiseData(df, dt, sym(variable))
  }
  
  dist_p <- dt %>%
    dplyr::filter(!grepl("Missing", !!sym(pair_variable))) %>%
    ggplot(aes(x=1-rmNA, color=!!sym(pair_variable))) +
    ggdist::stat_slab(alpha = .3) +
    stat_pointinterval(position = position_dodge(width = .4, preserve = "single")) +
    labs(
      x="Mitochondrial genetic similarity\n(Excluding missing positions)",
      y = NULL,
      color = variable_title) +
    scale_y_continuous(breaks = NULL)
  if(!is.null(color_vector)){
    dist_p <- dist_p + scale_color_manual(values = color_vector)
    }
  
  joint_p <- ggpubr::ggarrange(count_p, dist_p, nrow=1, align = "h", widths = c(1,1.5))
  return(joint_p)
} 

################################################################################
# multipanel plots
ClusterPlots <- function(df, cluster_id=NA, output_name=NA, filter_min=0, pairwise_df = diff_all, 
                         # sample_name="wormnum_dpdxid"
                         sample_name = "sample"){
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
    host <- df_cluster %>% ggplot(aes(x=as.character(year), y=after_stat(count), fill=host))+
      geom_bar(colour="black") +
      scale_fill_manual(values = host_colors, name="Host") +
      labs(x="Year", y="Specimens") + 
      theme(legend.position="top") +
      ggtitle("")
    barcodes <- df_cluster %>% ggplot(aes(x=as.character(year), y=after_stat(count), fill=amplicon))+
      geom_bar(colour="black") +
      scale_fill_manual(values = nextstrain_colors,
                        name="Barcode\ngroup") +
      labs(x="Year", y="Specimens") + 
      theme(legend.position="top") +
      guides(fill=guide_legend(nrow=2, byrow=TRUE))
    
    general_joint <- ggpubr::ggarrange(host, barcodes, nrow = 2, widths = c(1,1), align="hv")
  
    # relatedness and missingness plot
    # check the new sample_id is in the pairwise df
    if(!all(c(paste0(sample_name, ".x"), paste0(sample_name, ".y")) %in% names(pairwise_df))){
      print(paste("Adding", sample_name, "column to pairwise dataframe."))
      lookup_table <- setNames(df[[sample_name]], df[["sample"]])
      pairwise_df <- pairwise_df[, (paste0(sample_name, ".x")) := lookup_table[sample.x]]
      pairwise_df <- pairwise_df[, ("vassar_worm.x") := gsub(".*._", "", sample.x)]
      pairwise_df <- pairwise_df[, (paste0(sample_name, ".y")) := lookup_table[sample.y]]
      pairwise_df <- pairwise_df[, ("vassar_worm.y") := gsub(".*._", "", sample.y)]
    }

    tmp_pair <- pairwise_df %>% mutate(outline = ifelse(rmNA == 0, T, NA))
    related_joint <- related_and_missing_plot(tmp_pair, eval(sample_name), filter_min=filter_min)
    
    # all combined
    plot_name <- paste0("Cluster ID: ", output_name, "\n", 
                        "Narrative samples: ", nrow(df_cluster), "; Sequenced specimens: ", sum(!is.na(df_cluster$amplicon)))
    plot <-  #ggpubr::ggarrange(
      ggpubr::ggarrange(general_joint, related_joint, ncol=2, widths = c(1,2))
    # position, nrow=2, heights = c(2.75,1))
    plot <- ggpubr::annotate_figure(plot, top = plot_name)

    # save plot
    plot_height <- ifelse(nrow(df_cluster) <= 15, 5, 7.5) 
    plot_width <- ifelse(nrow(df_cluster) <= 15, 8, 12)
    SavePlots(plot, output_dir, paste0(output_name, ".png"), height=plot_height, width=plot_width)
    
    potential_outliers <- "Exit"
  } else{
    potential_outliers <- "Skipped"
  }
  return(plot)
} 

################################################################################
# Saving/exporting
################################################################################
SavePlots <- function(ggplot_obj, plot_path, name, width=8, height=4, dpi=100){
  ggsave(paste0(format(Sys.time(), "%Y%m%d"), "_", name), 
         plot = ggplot_obj,
         path = plot_path,
         width = width, height = height, units = c("in"), 
         dpi = dpi)
}