library(ggplot2)
library(tidyverse)

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "input.csv", help = "Input CSV file path"),
  make_option(c("-p", "--plot"), type = "character", default = NA, help = "Output plot file path. Ending is used to give file type as supported by ggplot2. (Ex: pdf, png)"),
  make_option(c("-o", "--output"), type = "character", default = NA, help = "Output CSV file path of top 20 positions with the most soft-clipping.")
)

opt_parser <- OptionParser(usage = "Usage: %prog [options]", option_list = option_list)

# Parse the command-line arguments
opt <- parse_args(opt_parser)

input_path <- opt$input
plot_path <- opt$plot
output_path <- opt$output

X = read_csv(input_path)

# Output empty files
if (nrow(X) == 0) 
{
  Y <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(Y) = c("n", "bin_coord", "bin_start", "bin_end", "direction")

  p = ggplot(Y, aes(x=bin_coord, y=n, color=direction)) +
  geom_point() +
  theme_bw()
  ggsave(plot_path, plot=p, width=12, height=6)

  Y = Y %>% ungroup() %>% select(n, bin_start, bin_end, direction)
  write_csv(Y, output_path)

  quit() 
}

bin_size = 10
X$bin = floor((X$position-1)/bin_size)+1
Y = X %>% group_by(bin, direction) %>% summarize(n=n())

Y$bin_coord = (Y$bin - 1) * bin_size
Y$direction = as.factor(Y$direction)

Y = Y %>% filter(n>1)

if (!is.na(plot_path)) {
  p = ggplot(Y, aes(x=bin_coord, y=n, color=direction)) +
    geom_point() +
    theme_bw()
  
  ggsave(plot_path, plot=p, width=12, height=6)
}


if (!is.na(output_path)) {
  Y = Y %>% arrange(desc(n), bin_coord) %>% head(20)
  Y$bin_start = (Y$bin - 1) * bin_size + 1
  Y$bin_end = (Y$bin) * bin_size
  Y = Y %>% ungroup() %>% select(n, bin_start, bin_end, direction)
  write_csv(Y, output_path)
}

