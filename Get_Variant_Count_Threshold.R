library(tidyverse)
library(data.table)

# args
# 1 - concatenated population count file
# 2 - percent threshold fr classifying outliers
# example
# args <- c("~/transfer/test_select.tsv", "99")

args <- commandArgs(TRUE)

all_variants <- data.table::fread(args[1], header = T) %>%
  dplyr::select(CHROM, START_BIN, STRAIN, COUNT)

df_dist_COUNT <- all_variants  %>%
  dplyr::arrange(COUNT) %>%
  dplyr::mutate(rownum = as.numeric(rownames(.))) %>%
  dplyr::group_by(COUNT) %>%
  dplyr::mutate(rank_COUNT=min(rownum)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rank_percent = ifelse(COUNT == 0, 0, rank_COUNT/nrow(all_variants) * 100))

df_dist_COUNT_99 <- df_dist_COUNT %>%
  dplyr::distinct(COUNT, .keep_all = T) %>%
  dplyr::filter(rank_percent >= as.numeric(args[2]))

COUNT_threshold <- min(df_dist_COUNT_99$COUNT)

write.table(COUNT_threshold, file = "Count_Threshold.tsv", col.names = F, row.names = F, sep = "\t")
