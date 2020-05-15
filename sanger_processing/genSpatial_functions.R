# tf_colors <- setNames(wesanderson::wes_palette("Cavalcanti1")[4:5], c("TRUE", "FALSE"))
tf_colors <- setNames(c("black", "grey50"), c("True", "False"))

distDensity_plots <- function(df, plot_name, project_dir){
  cdf <- ggplot(df, aes(x=dist_m/1000)) + 
    stat_ecdf(aes(colour = bc_match), size=2) +
    labs(x="Distance (km)", y="Cumulative density") +
    scale_color_manual(values=tf_colors, name="Same\nbarcode")
  pdf <-  ggplot(df, aes(x=dist_m/1000)) + 
    stat_density(aes(fill = bc_match), position = "identity", alpha=0.75) +
    scale_fill_manual(values=tf_colors, name="Same\nbarcode") +
    labs(x="Distance (km)", y="Probability density")
  #df_plots <- ggpubr::ggarrange(pdf, cdf, common.legend=TRUE, legend = "top")
  df_plots <- ggpubr::ggarrange(pdf, #+ theme(axis.title=element_blank(), axis.text=element_blank(), plot.margin = unit(c(2,2,2,2), "lines")), 
                                cdf, #+ theme(axis.title=element_blank(), axis.text=element_blank(), plot.margin = unit(c(2,2,2,2), "lines")), 
                                common.legend=TRUE, legend = "top")
  ggsave(paste0(plot_name, ".png"), plot = df_plots, path = paste(project_dir, "plots", sep="/"), 
         width = 7, height = 3.5, units = c("in"), dpi = 300)
  return(df_plots)
}
  

# modified from https://www.r-bloggers.com/extract-different-characters-between-two-strings-of-equal-length/
string_diff <- function(a, b){
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  seq.a <- unlist(strsplit(a,split=""))
  seq.b <- unlist(strsplit(b,split=""))
  diff.d <- rbind(seq.a,seq.b)
  only.diff <-diff.d[,diff.d[1,]!=diff.d[2,]]
  pos <- which(diff.d[1,]!=diff.d[2,])
  if(0 < length(pos) & length(pos) <= 2){
    return(cbind(full_barcode.x=a, full_barcode.y=b, diff_pos=toString(pos)))
  }else{
    return(cbind(full_barcode.x=a, full_barcode.y=b, diff_pos="ignore"))
  }
}