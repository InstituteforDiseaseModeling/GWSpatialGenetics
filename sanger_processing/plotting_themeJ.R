theme_j <- function () {
  theme_bw(base_size=16) %+replace%
    theme(
      # font sizes and color
      panel.background  = element_blank(),
      plot.background   = element_rect(fill="transparent", colour=NA),
      plot.title        = element_text(size = rel(.85)),
      strip.background  = element_rect(fill="transparent", colour=NA),
      strip.text        = element_text(face="bold", size=rel(.6)),
      axis.title        = element_text(size=rel(0.8)),
      axis.text         = element_text(size=rel(0.6), color="grey30"),
      # legend
      legend.title         = element_text(size=rel(0.8)),
      legend.text          = element_text(size=rel(0.6)),
      legend.background    = element_rect(fill="transparent", colour=NA),
      legend.key           = element_rect(fill="transparent", colour=NA),
      legend.justification = "top"
    )
}

theme_set(theme_j())