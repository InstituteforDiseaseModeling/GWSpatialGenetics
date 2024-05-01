################################################################################
# Purpose: Sankey plots for different sequencing technologies
# Author: Jessica Ribado
# Date: 01/2022
################################################################################

################################################################################
# set-up 
################################################################################
# load libraries
for(p in c('vcfR', 'data.table', 'dplyr', 'tidyr', "ggplot2",
           'ggsci', 'ggforce', 'ggpubr', 'knitr')){
  if(!p %in% installed.packages()[,1]){
    print(p)
    install.packages(p, repos =  "https://cloud.r-project.org", dependencies = T )
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

# set global options
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE)


# functions
sankey_plot <- function(df){
  
  # set order of barcodes for plotting
  otu_ord <- gtools::mixedsort(unique(df$amplicon))
  # otu_color <- paste0('d3.scaleOrdinal() .domain([', 
  #                     paste0('"', paste(ggsci_cols$cluster, collapse='", "'), '"'),']) .range([', 
  #                     paste0('"', paste(ggsci_cols$color, collapse='", "'), '"'), '])')
  
  # plot
  nodes <- data.frame(
    new_barcode=c(as.character(df$amplicon), 
                  as.character(df$kinship_group)) %>% unique())
  
  df$IDsource <- match(df$amplicon, nodes$new_barcode)-1
  df$IDtarget <- match(df$kinship_group, nodes$new_barcode)-1
  
  # Make the Network
  p <- sankeyNetwork(Links = df, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "count", NodeID = "new_barcode", 
                     #colourScale=noquote(otu_color),
                     sinksRight=FALSE)
  p
  
  file = paste0(format(Sys.time(), "%Y%m%d"), "_", "202108_AmpliconKinshipGroups.png")
  setwd(out_path)
  saveNetwork(p, file=file, selfcontained = TRUE)
  webshot(file, gsub("html", "png", file), vwidth = 1000, vheight = 900)
}

################################################################################
# load and format data
################################################################################
out_path = '/mnt/data/guinea_worm/analyses'
seq_groups <- read.delim('/mnt/data/guinea_worm/nextstrain/round_202108_known/metadata_clusters.tsv')

df <- dplyr::filter(seq_groups, !is.na(kinship)) %>%
  group_by(amplicon, kinship_group) %>%
  summarise(count = n())

sankey_plot(df)
