#' produces a figure summarizing alpha and beta diversity on scatterplots
#'
#' @param alpha numeric. alpha diversity estimated
#' @param PCoA_1 numeric. first component obtained from PCoA
#' @param PCoA_2 numeric. second component obtained from PCoA
#' @param PCoA_3 numeric. third component obtained from PCoA
#' @param classes character. classes of vegetation types/species communities
#' @param filename character. path where to save figure
#' @param palette character. palette to be used.
#' --> either name for a palette in RColorBrewer::brewer.pal or set of colors
#'
#' @return none
#' @export

scatter_alphabeta <- function(alpha, PCoA_1, PCoA_2, PCoA_3, classes,
                              filename = NULL, palette = 'Set1'){

  # define colors corresponding to different classes
  nbClasses <- length(unique(classes))
  if (length(palette) == 1) cols <- RColorBrewer::brewer.pal(n = nbClasses, name = palette)
  if (length(palette) > 1 & length(palette) == nbClasses) {
    cols <- palette
  } else if (length(palette) > 1 & !length(palette) == nbClasses) {
    message('warning in function "scatter_alphabeta"')
    message('please specify a number of colors in palette corresponding to the number of classes')
    stop()
  }
  # produce a data.frame corresponding to variables to be used
  diversity_summary <- data.frame('alpha' = alpha,
                                  'PCoA_1' = PCoA_1,
                                  'PCoA_2' = PCoA_2,
                                  'PCoA_3' = PCoA_3,
                                  'classes' = classes)
  # produce subplots coresponding to pairwise comparison of PCoA components
  aes <- theme <- NULL
  g1 <- ggplot2::ggplot(diversity_summary,
                        aes (x = PCoA_1, y = PCoA_2, color = classes, size = alpha)) +
    ggplot2::geom_point(alpha=0.6) +
    ggplot2::scale_color_manual(values=cols) + theme(legend.position="bottom")
  g2 <- ggplot2::ggplot(diversity_summary,
                        aes (x = PCoA_1, y = PCoA_3, color = classes, size = alpha)) +
    ggplot2::geom_point(alpha=0.6) +
    ggplot2::scale_color_manual(values=cols)
  g3 <- ggplot2::ggplot(diversity_summary,
                        aes (x = PCoA_2, y = PCoA_3, color = classes, size = alpha)) +
    ggplot2::geom_point(alpha=0.6) +
    ggplot2::scale_color_manual(values=cols)
  #extract legend
  #https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  get_legend <- function(a.gplot){
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  # produce full plot arranging subplots
  legend <- get_legend(g1)
  gAll <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(g1 + theme(legend.position="none"),
                           g2 + theme(legend.position="none"),
                           g3 + theme(legend.position="none"),
                           nrow=1), legend, nrow=2, heights=c(8, 1))

  if (!is.null(filename))
    ggplot2::ggsave(filename, plot = gAll, device = 'png', scale = 1, width = 15,
           height = 5, units = "in", dpi = 600, limitsize = TRUE)
  return(invisible())
}
