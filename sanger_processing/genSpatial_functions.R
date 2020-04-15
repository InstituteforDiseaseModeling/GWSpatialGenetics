tf_colors <- setNames(wesanderson::wes_palette("Cavalcanti1")[4:5], c("TRUE", "FALSE"))

distDensity_plots <- function(df, plot_name, project_dir){
  cdf <- ggplot(df, aes(x=dist_m/1000)) + 
    stat_ecdf(aes(colour = bc_match), size=2) +
    labs(x="Distance (km)", y="Cumulative density") +
    scale_color_manual(values=tf_colors, name="Same\nbarcode")
  pdf <-  ggplot(df, aes(x=dist_m/1000)) + 
    stat_density(aes(fill = bc_match), position = "identity", alpha=0.75) +
    scale_fill_manual(values=tf_colors, name="Same\nbarcode") +
    labs(x="Distance (km)", y="Probability density")
  df_plots <- ggpubr::ggarrange(cdf, pdf, common.legend=TRUE, legend = "top")
  ggsave(paste0(plot_name, ".png"), plot = df_plots, path = paste(project_dir, "plots", sep="/"), 
         width = 8, height = 4, units = c("in"), dpi = 300)
  return(df_plots)
}
  
