library(ggplot2)
library(RColorBrewer)

my.theme <- theme_light() + theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  plot.caption = element_text(size = 24),
  axis.title = element_text(size = 22),
  axis.text = element_text(size = 18),
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 14),
  strip.background = element_blank(),
  strip.text = element_text(size = 28, colour = "black")
)

#' Plot PCA data in a legible scatterplot, with experimentals
#' overlaid on top of reference 1KG subjects.
#'
#' @param in.ancestry.fn character; input somalier ancestry tsv file
#' @param out.plot.fn character; output PCA plot filename. expected png,
#' but with support for anything ggplot2/R can handle with ggsave.
#' @param plot.count.table logical; whether to include a count table
#' rendered in the plot itself.
#' @return NULL
run.plot.pca <- function(in.ancestry.fn, out.plot.fn, plot.count.table) {
  stopifnot(length(in.ancestry.fn) == 1, is.character(in.ancestry.fn))
  stopifnot(length(out.plot.fn) == 1, is.character(out.plot.fn))
  stopifnot(length(plot.count.table) == 1, is.logical(plot.count.table))

  df <- read.table(in.ancestry.fn,
    header = TRUE, comment.char = "",
    quote = "", sep = "\t", stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ## annotate input data with whether they are experimentals.
  df$data.source <- ifelse(df[, "given_ancestry"] == "", "Current study", "1000 Genomes")

  ## set factor level order
  df$data.source <- factor(df$data.source, levels = c("Current study", "1000 Genomes"))
  ## make a scatterplot with substantial formatting
  my.plot <- ggplot(aes(
    x = PC1, y = PC2, colour = predicted_ancestry,
    shape = data.source,
    alpha = data.source,
    size = data.source
  ), data = df)
  my.plot <- my.plot + geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3)
  my.plot <- my.plot + my.theme + geom_point()
  ## this palette is selected for contrast within a colorblind-friendly set.
  ## see other options with display.brewer.all(colorblindFriendly = TRUE).
  ## count assumes the five 1KG supercontinents.
  my.colors <- brewer.pal(6, "Dark2")[c(1:4, 6)]
  my.plot <- my.plot + scale_colour_manual(name = "Supercontinent", values = my.colors)
  my.plot <- my.plot + scale_shape_manual(name = "Data Source", values = c(15, 1))
  my.plot <- my.plot + scale_alpha_manual(values = c(1, 0.4))
  my.plot <- my.plot + scale_size_manual(values = c(2, 1))
  my.plot <- my.plot + guides(alpha = "none", size = "none")
  if (plot.count.table) {
    my.plot <- my.plot + annotate("segment", x = 50.5, xend = 59, y = 38.25, yend = 38.25, alpha = 0.5)
    my.plot <- my.plot + annotate("segment", x = 55, xend = 55, y = 30.75, yend = 39.75, alpha = 0.5)
    my.plot <- my.plot + annotate("text", x = c(52.5, 57.5), y = 39, label = c("Ancestry", "Count"))
    anc.counts <- table(as.character(df[df$data.source == "Current study", "predicted_ancestry"]))
    my.plot <- my.plot + annotate("text", x = 52.5, y = seq(37.5, 31.5, length.out = 5), label = names(anc.counts))
    my.plot <- my.plot + annotate("text", x = 57.5, y = seq(37.5, 31.5, length.out = 5), label = unname(anc.counts))
  }
  ggsave(out.plot.fn, plot = my.plot, height = 10, width = 16 / 9 * 10, units = "in")
}

if (exists("snakemake")) {
  run.plot.pca(
    snakemake@input[["somalier_ancestry"]],
    snakemake@output[["plotname"]],
    FALSE
  )
}
