library(tidyverse)
library(smoother)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

outliers <- readr::read_tsv("CB4856_Outlier_Counts.tsv")
coverage <- readr::read_tsv("CB4856.regions.bed.gz", col_names = c("CHROM", "START_BIN", "END_BIN", "COVERAGE"))

# sv file - DEL
# awk '$5 == "DEL" {print}' CB4856_HQ_SnpEff.bed | cut -f -3 | bedtools coverage -a 1000windows.bed -b stdin > CB4856_DEL_nSVs.bed
# awk '$5 == "INV" {print}' CB4856_HQ_SnpEff.bed | cut -f -3 | bedtools coverage -a 1000windows.bed -b stdin > CB4856_INV_nSVs.bed
# sv file
# grep -v BND CB4856_HQ_SnpEff.bed | cut -f -3 | bedtools coverage -a 1000windows.bed -b stdin > CB4856_nSVs.bed

# 
# grep -v BND CB4856_HQ_SnpEff.bed | cut -f -3 | awk '$3-$2 < 1e5 {print}' | bedtools coverage -a 1000windows.bed -b stdin > CB4856_nSVs.bed
# 
# 
# grep -v BND CB4856_SnpEff.bed | cut -f -3 | awk '$3-$2 < 1e5 {print}' | bedtools coverage -a 1000windows.bed -b stdin > CB4856_nSVs_low.bed

total_svs <- readr::read_tsv("CB4856_nSVs_low.bed", col_names = c("CHROM", "START_BIN", "END_BIN", "total_SV_CT","total_SV_bases","bin_size","fraction_SV_bases")) %>%
  dplyr::select(-bin_size)

smooth_outs <- outliers %>%
  dplyr::left_join(., coverage,  by = c("CHROM", "START_BIN", "END_BIN")) %>%
  dplyr::left_join(., total_svs,  by = c("CHROM", "START_BIN", "END_BIN")) %>%
  dplyr::filter(CHROM=="V")%>%
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
  dplyr::ungroup() %>%
  dplyr::mutate(med_cov = median(COVERAGE),
                iqr_cov = IQR(COVERAGE)) %>%
  dplyr::mutate(iqr_cov_out = ifelse(COVERAGE < med_cov-(iqr_cov), "low", 
                                     ifelse(COVERAGE > med_cov+(iqr_cov), "high", "not_outlier"))) %>%
  dplyr::mutate(zero_cov = ifelse(iqr_cov_out == "low" , 100, 1)) %>%
  dplyr::mutate(region = paste0(CHROM, ":", START_BIN, "-", END_BIN)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(smoothed_cov = smth.gaussian(zero_cov, 
                                              alpha = 5,
                                              window = 10,
                                              tails = T)) %>%
  dplyr::mutate(norm_smooth_out = smoothed_outs/max(smoothed_outs),
                norm_coverage = log(smoothed_cov)/max(log(smoothed_cov)),
                norm_svct = total_SV_CT/max(total_SV_CT))

ggplot()+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  # geom_line(aes(y=COVERAGE,x = MID_BIN/1e6), color = "purple", data = smooth_outs_pr) +
  geom_line(aes(y=norm_smooth_out,x = MID_BIN/1e6), color = "red", data = smooth_outs_pr) +
  geom_line(aes(y=-norm_coverage,x = MID_BIN/1e6), color = "cyan", data = smooth_outs_pr) +
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank()) +
  labs(color = "Outlier\nRank") + xlim(c(15,20))


make_windows <- smooth_outs_pr %>%
  dplyr::select(CHROM:MID_BIN, COVERAGE, norm_smooth_out, norm_coverage, norm_svct) %>%
  dplyr::mutate(norm_coverage = ifelse(norm_coverage < 0, 0, norm_coverage))

ggplot()+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  geom_line(aes(y=smoothed_outs,x = MID_BIN/1e6), color = "red", data = smooth_outs_pr) +
  # geom_line(aes(y=COVERAGE,x = MID_BIN/1e6), color = "cyan", data = smooth_outs_pr) +
  geom_line(aes(y=total_SV_CT,x = MID_BIN/1e6), color = "purple", data = smooth_outs_pr) +
  theme_bw(18) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank()) +
  labs(color = "Outlier\nRank") + xlim(2,3)


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
