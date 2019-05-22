#!/usr/bin/env Rscript
library(tidyverse)
library(anomalize)
library(tibbletime)

# args 
# 1 - sample name
# 2 - count outlier df
# 3 - coverage bed file
# 4 - sv count bed file
args <- commandArgs(trailingOnly = TRUE)

outliers <- readr::read_tsv(args[2])
coverage <- readr::read_tsv(args[3], col_names = c("CHROM", "START_BIN", "END_BIN", "COVERAGE"))

total_svs <- readr::read_tsv(args[4], col_names = c("CHROM", "START_BIN", "END_BIN", "total_SV_CT","total_SV_bases","bin_size","fraction_SV_bases")) %>%
  dplyr::select(-bin_size)

region_info <- outliers %>%
  dplyr::left_join(., coverage,  by = c("CHROM", "START_BIN", "END_BIN")) %>%
  dplyr::left_join(., total_svs,  by = c("CHROM", "START_BIN", "END_BIN")) 

for(chrom in unique(region_info$CHROM)){
  
  smooth_outs <- region_info %>%
    dplyr::filter(CHROM==chrom)%>%
    dplyr::mutate(direction = ifelse(is.na(direction), "ignore_na", 
                                     ifelse(direction == "Down", "ignore_down", "Up"))) %>%
    dplyr::mutate(zero_counts = ifelse(outlier == "Yes", COUNT, 0)) %>%
    dplyr::mutate(smoothed_outs = smth.gaussian(zero_counts, 
                                                alpha = 5,
                                                window = 10,
                                                tails = T))
  # filter wholegenome dataset to chrom
  cov_oulier <- smooth_outs %>%
    dplyr::select(CHROM:MID_BIN,COVERAGE) %>%
    dplyr::ungroup()

  # make fake date to use tibbletime and anomalize functionality
  fake_time <- tibbletime::create_series('1900' ~ '2017', 'daily', class = "Date")[1:nrow(cov_oulier),]
  
  # find outliers
  temp_outlier_time <- cov_oulier %>%
    dplyr::mutate(date = fake_time$date,
                  COVERAGE = ifelse(COVERAGE == 0, .01, log(COVERAGE))) %>%
    as.tbl() %>%
    time_decompose(COVERAGE, method = "stl", trend = 50) %>%
    anomalize(remainder, method = "gesd", alpha = 0.005, max_anoms = 0.05, verbose = T)
  
  # extract outlier information
  outlier_df <- temp_outlier_time$anomaly_details$outlier_report %>%
    dplyr::select(index, cov_rank = rank, cov_outlier = outlier, cov_direction = direction)
  
  # append outlier information count df
  outlier_by_chrom <- smooth_outs %>%
    dplyr::ungroup() %>%
    dplyr::mutate(index = 1:n()) %>%
    dplyr::left_join(.,outlier_df, by = "index")
  
  # deal with NAs from join
  outlier_by_chrom$cov_outlier[is.na(outlier_by_chrom$cov_outlier)] <- "No"
  outlier_by_chrom$cov_direction[is.na(outlier_by_chrom$cov_direction)] <- "NA"
  
  # append chromosomal region information
  if(!exists("count_regions")){
    count_regions <- outlier_by_chrom %>%
      dplyr::select(CHROM:END_BIN, MID_BIN, index, GENOMIC_REGION, STRAIN, COUNT, count_rank = rank, count_outlier = outlier, count_direction = direction, smoothed_count_outlier = smoothed_outs,
                    COVERAGE, cov_rank:cov_direction, total_SV_CT:fraction_SV_bases)
  } else {
    temp_regions <- outlier_by_chrom %>%
      dplyr::select(CHROM:END_BIN, MID_BIN, index, GENOMIC_REGION, STRAIN, COUNT, count_rank = rank, count_outlier = outlier, count_direction = direction, smoothed_count_outlier = smoothed_outs,
                    COVERAGE, cov_rank:cov_direction, total_SV_CT:fraction_SV_bases)
    count_regions <- dplyr::bind_rows(count_regions, temp_regions)
  }
}
# make non-outlier bins rank to 1+the max outlier rank (for plotting purposes)
count_regions$cov_rank[is.na(count_regions$cov_rank)] <- max(count_regions$cov_rank,na.rm = T)+1


write.table(count_regions, 
            file = glue::glue("{args[1]}_Processed_Outliers.tsv"), 
            quote = F, col.names = T, row.names = F, sep = "\t")
