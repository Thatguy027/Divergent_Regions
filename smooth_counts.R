library(tidyverse)
library(smoother)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

outliers <- readr::read_tsv("CB4856_Outlier_Counts.tsv")
coverage <- readr::read_tsv("CB4856.regions.bed.gz", col_names = c("CHROM", "START_BIN", "END_BIN", "COVERAGE"))

smooth_outs <- outliers %>%
  dplyr::left_join(., coverage,  by = c("CHROM", "START_BIN", "END_BIN")) %>%
  dplyr::filter(CHROM=="II")%>%
  dplyr::mutate(direction = ifelse(is.na(direction), "ignore_na", 
                                   ifelse(direction == "Down", "ignore_down", "Up"))) %>%
  dplyr::mutate(zero_counts = ifelse(outlier == "Yes", COUNT, 0)) %>%
  dplyr::mutate(smoothed_outs = smth.gaussian(zero_counts, 
                                              alpha = 5,
                                              window = 10,
                                              tails = T))

smooth_outs_pr <- smooth_outs %>%
  dplyr::group_by(GENOMIC_REGION) %>%
  dplyr::mutate(start_region = min(START_BIN),
                end_region = max(END_BIN),
                med_count = mean(COUNT),
                iqr_count = sd(COUNT)) %>%
  dplyr::mutate(count_out = ifelse(COUNT > med_count+(7*iqr_count), "7iqr",
                                 ifelse(COUNT > med_count+(6*iqr_count), "6iqr",
                                        ifelse(COUNT > med_count+(5*iqr_count), "5iqr",
                                               ifelse(COUNT > med_count+(4*iqr_count), "4iqr",
                                                      ifelse(COUNT > med_count+(3*iqr_count), "3iqr",
                                                             ifelse(COUNT > med_count+(2*iqr_count), "2iqr",
                                                                    ifelse(COUNT > med_count+(1*iqr_count), "1iqr", "not_outlier")))))))) %>%
  dplyr::select(-index, -med_count, -iqr_count) %>%
  dplyr::mutate(med_cov = median(COVERAGE),
                iqr_cov = IQR(COVERAGE)) %>%
  dplyr::mutate(iqr_cov_out = ifelse(COVERAGE < med_cov-(iqr_cov), "low", 
                                     ifelse(COVERAGE > med_cov+(iqr_cov), "high", "not_outlier"))) %>%
  dplyr::mutate(region = paste0(CHROM, ":", START_BIN, "-", END_BIN))

ggplot()+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  # geom_point( aes(color = rank, y = zero_counts,x = MID_BIN/1e6), data = smooth_outs) +
  geom_line(aes(y=smoothed_outs,x = MID_BIN/1e6), color = "black", data = smooth_outs_pr) +
  geom_segment(aes(x=start_region/1e6, xend = end_region/1e6, y = med_count+(3*iqr_count), yend = med_count+(3*iqr_count), color = GENOMIC_REGION), size =2,
               data = smooth_outs_pr %>% dplyr::distinct(GENOMIC_REGION, end_region, start_region, .keep_all=T)) +
  # geom_rect(aes(xmin=START_DIV/1e6, xmax=END_DIV/1e6, ymin=-5, ymax=0), data = div_chrom, fill = "cyan") +
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) +
  labs(color = "Outlier\nRank") 

ggplot()+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  scale_color_viridis_c(direction = -1, option = "B")+
  geom_point( aes(color = rank, y = COUNT,x = MID_BIN/1e6), data = smooth_outs) +
  geom_line(aes(y=smoothed_outs,x = MID_BIN/1e6), color = "red", data = smooth_outs) +
  # geom_rect(aes(xmin=START_DIV/1e6, xmax=END_DIV/1e6, ymin=-5, ymax=0), data = div_chrom, fill = "cyan") +
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) +
  labs(color = "Outlier\nRank") 

smooth_outs %>%
  dplyr::mutate( masked = ifelse(smoothed_outs > 0 & COUNT > 150, "MASKED", "PASS")) %>%
  ggplot()+
  aes(x = COVERAGE, fill = masked) +
  geom_density(alpha = 0.5)+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) +
  labs(fill = "Region") +
  xlim(c(0, 100))

ggplot(smooth_outs)+
  aes(x = smoothed_outs, y = COVERAGE, color = rank) +
  geom_point()+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  scale_color_viridis_c(direction = -1, option = "B")+
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Variant Count", y = "Coverage")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) +
  labs(color = "Outlier\nRank") +
  ylim(c(0, 100))

smooth_outs %>%
  dplyr::filter(COVERAGE < 600) %>%
  dplyr::mutate(norm_cov = scale(COVERAGE)) %>%
  dplyr::mutate(norm_out = scale(smoothed_outs)) %>%
  ggplot()+
  facet_grid(CHROM~., space = "free", scales = "free")+
  scale_color_viridis_c(direction = -1, option = "B")+
  geom_line(aes(y=norm_cov,x = MID_BIN/1e6), color = "cyan") +
  geom_line(aes(y=norm_out,x = MID_BIN/1e6), color = "red") +
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) +
  labs(color = "Outlier\nRank") +
  ylim(-3,3)

smooth_outs %>%
  dplyr::mutate(masked = ifelse(smoothed_outs > 0, "masked", "pass")) %>%
  dplyr::filter(masked == "masked") %>%
  dplyr::mutate(masked_size = END_BIN-START_BIN) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(masked_region = sum(masked_size))

smooth_outs %>%
  dplyr::mutate(masked = ifelse(smoothed_outs > 0, "masked", "pass")) %>%
  dplyr::filter(masked == "masked") %>%
  dplyr::mutate(masked_size = END_BIN-START_BIN) %>%
  # dplyr::group_by(CHROM) %>%
  dplyr::summarise(masked_region = sum(masked_size))

smooth_outs %>%
  dplyr::mutate( masked = ifelse(smoothed_outs > 0 & COUNT > 150, "MASKED", "PASS")) %>%
  dplyr::group_by(masked) %>%
  dplyr::summarise(masked_region = sum(COUNT))
