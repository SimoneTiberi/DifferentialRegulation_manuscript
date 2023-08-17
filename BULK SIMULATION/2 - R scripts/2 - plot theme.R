library(RColorBrewer);
all_colours = c(
  "#F4A582", # "DEXSeq_ECs"
  "#A63603",  # "DEXSeq_TECs"
  "#1F78B4", # "DifferentialRegulation" #CAB2D6
  "#66A61E", # "DRIMSeq"
  "#A6761D",  # eisaR
  "#404040", # "satuRn"
  "#FFD92F" # "SUPPA2"
)

methods_all = c(
  "DEXSeq_ECs",
  "DEXSeq_TECs",
  "DifferentialRegulation",
  "DRIMSeq",
  "eisaR",
  "satuRn",
  "SUPPA2"
)

colors_method = all_colours
names(colors_method) = methods_all

# points, borders:
shape_border = c("DEXSeq_ECs" = 0,
                 "DEXSeq_TECs" = 1,
                 "DifferentialRegulation" = 2,
                 "DRIMSeq" = 3,
                 "eisaR" = 4,
                 "satuRn" = 5,
                 "SUPPA2" = 6)

# points, fill:
shape_fill = c("DEXSeq_ECs" = 3,
               "DEXSeq_TECs" = 15,
               "DifferentialRegulation" = 16,
               "DRIMSeq" = 17,
               "eisaR" = 23,
               "satuRn" = 25,
               "SUPPA2" = 26)

.prettify <- function(theme = NULL, ...) {
  if (is.null(theme)) theme <- "classic"
  base <- paste0("theme_", theme)
  base <- getFromNamespace(base, "ggplot2")
  base(base_size = 8) + theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(2, "mm"),
    strip.background = element_rect(fill = NA),
    plot.margin = unit(rep(1, 4), "mm"),
    panel.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,1,0,"mm"),
    ...)}

my_theme <- theme(strip.background = element_blank(),
                  strip.text = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(4)),
                  axis.text.y=element_text(size=rel(4)),
                  axis.title.y = element_text(size=rel(4)),
                  axis.title.x = element_text(size=rel(4)),
                  legend.title=element_text(size=rel(1.5)),
                  legend.text=element_text(size=rel(4)),
                  legend.key.width=unit(2, "cm"),
                  aspect.ratio = 1, legend.position="bottom",
                  legend.box="vertical", legend.margin=margin())
