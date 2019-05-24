library(tidyverse)
library(zoo)
library(rlang)

options(scipen=999)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# args
# 1 - processed outlier file
# 2 - window size (--masked_range)
# 3 - outlier rank for masking (--outlier_rank)
# 4 - variant count for masking (--outlier_count) # daehan determins this
# 5 - fraction of window to be overlapped by SV (--sv_fraction)
# 6 - minimum variant count to be considered diverged (--min_outlier_count)
# example
# args <- c("CB4856_Processed_Outliers.tsv","10","100","50", "0.5", "20")

args <- commandArgs(TRUE)

leads <- function(var, n=as.numeric(args[2])){
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(dplyr::lead(!!var, !!.x)) ) %>% 
    set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
  
}

lags <- function(var, n=as.numeric(args[2])){
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(dplyr::lag(!!var, !!.x)) ) %>% 
    set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
  
}


outliers <- readr::read_tsv(args[1])




# Coverage - remove coverage outliers

m_windows <- outliers %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, .keep_all = T) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(window_mask = ifelse(cov_direction %in% c("Down"), "Masked_Low_Coverage", "Pass")) %>% # coverage mask
  dplyr::mutate(window_mask = ifelse(cov_direction %in% c("Up"), "Masked_High_Coverage", window_mask)) %>% # coverage mask
  dplyr::mutate(window_mask = ifelse(count_rank < as.numeric(args[3]), "Masked_Outlier", window_mask)) %>% # outlier mask
  # dplyr::mutate(window_mask = ifelse(COUNT > as.numeric(args[4]), "Masked_Count", window_mask)) %>% # outlier mask
  dplyr::mutate(window_mask = ifelse(count_outlier == "Yes" & fraction_SV_bases > as.numeric(args[5]), "Masked_SV", window_mask)) %>% # High SV fraction and variant count outlier mask
  # If regions flanking are both masked, then mask central region
  dplyr::mutate(window_mask = ifelse((grepl("Masked", lag(window_mask)) & grepl("Masked", lead(window_mask))) & !grepl("Masked", window_mask), "Masked_Two_Flank", window_mask)) %>%  
  # If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse((lag(window_mask) == "Masked" | lead(window_mask) == "Masked") & count_outlier == "Yes", "Masked_One_Flank_Outlier", window_mask)) %>%
  # Second round of masking - If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", lag(window_mask)) | grepl("Masked", lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P3",window_mask)) %>%
  # Third round of masking - If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", lag(window_mask)) | grepl("Masked", lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P4",window_mask)) %>%
  # Fourth round of masking - If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", lag(window_mask)) | grepl("Masked", lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P5",window_mask)) %>%
  # Exentend masked regions by 10 windows
  # First go in leading direction
  dplyr::mutate(!!!leads(window_mask)) %>%
  dplyr::ungroup() %>%
  # Make long
  tidyr::gather(lead_step, mask, -(CHROM:window_mask)) %>%
  # Convert column names from leads function to numeric number of windows away
  dplyr::mutate(lead_step = as.numeric(factor(gsub("lag_window","", lead_step)))) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  dplyr::group_by(CHROM,START_BIN,END_BIN) %>%
  dplyr::distinct(mask,.keep_all=T) %>%
  # If there is only one mask flag per group, keep it. If more than one - mark as leading mask (deals with duplicates that result from gather function)
  dplyr::mutate(lead_mask = ifelse(n()==1, mask,
                                   ifelse(n() > 1, "Masked_Lead"))) %>%
  # Take distinct chrom, start bin, while maintaining the distance to the original masked region
  dplyr::arrange(CHROM, START_BIN, desc(lead_step)) %>%
  dplyr::distinct(lead_mask, .keep_all = T) %>%
  dplyr::rowwise() %>%
  # make leading mask column, either original mask or "Lead_Mask"
  dplyr::mutate(lead_mask = ifelse(grepl("Mask", window_mask), window_mask, lead_mask)) %>%
  dplyr::select(-mask) %>%
  dplyr::group_by(CHROM) %>%
  # Repeat above operation for other side of window
  dplyr::mutate(., !!!lags(window_mask)) %>%
  dplyr::ungroup() %>%
  tidyr::gather(lag_step, mask, -(CHROM:lead_mask)) %>%
  dplyr::mutate(lag_step = as.numeric(factor(gsub("lag_window","", lag_step)))) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  dplyr::group_by(CHROM,START_BIN,END_BIN) %>%
  dplyr::distinct(mask,.keep_all=T) %>%
  dplyr::mutate(lag_mask = ifelse(n()==1, mask,
                                  ifelse(n() > 1, "Masked_Lag")))%>%
  dplyr::arrange(CHROM, START_BIN, desc(lag_step)) %>%
  dplyr::distinct(lag_mask, .keep_all = T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(lag_mask = ifelse(grepl("Mask", window_mask), window_mask, lag_mask)) %>%
  dplyr::select(-mask) %>%
  dplyr::mutate(final_mask = ifelse(grepl("Mask",window_mask), window_mask,
                                    ifelse(grepl("Mask",lead_mask), lead_mask,
                                           ifelse(grepl("Mask", lag_mask), lag_mask, "Pass"))))


# Counts
m_windows_pr <- m_windows %>%
  dplyr::mutate(final_mask_loose = ifelse((final_mask %in% c("Masked_Lead", "Masked_Lag") & COUNT < as.numeric(args[6])) | (final_mask %in% c("Masked_Outlier") & COUNT < as.numeric(args[6]))  , "Pass", final_mask))

masked_plot <- ggplot(m_windows_pr)+
  facet_grid(CHROM~., space = "free", scales = "free")+
  geom_point(aes(x = MID_BIN/1e6, y = COUNT, color = final_mask)) +
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -5, yend = -5, color = final_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", final_mask_loose)))+
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -10, yend = -10, color = final_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", final_mask)))+
  theme_bw(18) +
  scale_color_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))

strain <- unique(m_windows$STRAIN)

ggsave(masked_plot, filename = glue::glue("{strain}_masked_plot.pdf"), height = 10, width = 24)

write.table(m_windows_pr, file = glue::glue("{strain}_Mask_DF.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

loose_mask <- m_windows_pr %>%
  dplyr::filter(grepl("Mask", final_mask_loose)) %>%
  dplyr::select(CHROM, START_BIN, END_BIN) %>%
  tidyr::unite(temp1, START_BIN,END_BIN, sep = "-") %>%
  tidyr::unite(REGION, CHROM, temp1, sep = ":")

write.table(loose_mask, file = glue::glue("{strain}_Loose_Mask.tsv"), col.names = F, row.names = F, quote = F, sep = "\t")

conservative_mask <- m_windows_pr %>%
  dplyr::filter(grepl("Mask", final_mask)) %>%
  dplyr::select(CHROM, START_BIN, END_BIN) %>%
  tidyr::unite(temp1, START_BIN,END_BIN, sep = "-") %>%
  tidyr::unite(REGION, CHROM, temp1, sep = ":")

write.table(conservative_mask, file = glue::glue("{strain}_Conservative_Mask.tsv"), col.names = F, row.names = F, quote = F, sep = "\t")


m_windows <- outliers %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, .keep_all = T) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(window_mask = ifelse(cov_direction %in% c("Down"), "Masked_Low_Coverage", "Pass")) %>% # coverage mask
  dplyr::mutate(window_mask = ifelse(cov_direction %in% c("Up"), "Masked_High_Coverage", window_mask)) %>% # coverage mask
  dplyr::mutate(window_mask = ifelse(count_direction == "Up", "Masked_Outlier", window_mask)) %>% # outlier mask
  dplyr::mutate(window_mask = ifelse(COUNT > as.numeric(args[4]) & !grepl("Masked", window_mask), "Masked_Count", window_mask)) %>% # outlier mask
  dplyr::mutate(window_mask = ifelse(count_outlier == "Yes" & fraction_SV_bases > as.numeric(args[5]), "Masked_SV", window_mask)) %>% # High SV fraction and variant count outlier mask
  dplyr::mutate(window_mask = ifelse( total_SV_CT > 2 & COUNT > 10, "Masked_SV_count", window_mask)) %>% # High SV fraction and variant count outlier mask
  dplyr::mutate(window_mask = ifelse( total_SV_CT > 5 & COVERAGE < 5, "Masked_SV_count_cov", window_mask)) %>% # High SV fraction and variant count outlier mask
  # If regions flanking are both masked, then mask central region
  dplyr::mutate(window_mask = ifelse((grepl("Masked", dplyr::lag(window_mask)) & grepl("Masked", dplyr::lead(window_mask))) & !grepl("Masked", window_mask), "Masked_Two_Flank", window_mask)) %>%  
  # If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse((dplyr::lag(window_mask) == "Masked" | dplyr::lead(window_mask) == "Masked") & count_outlier == "Yes", "Masked_One_Flank_Outlier", window_mask)) %>%
  # Second round of masking - If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", dplyr::lag(window_mask)) | grepl("Masked", dplyr::lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P3",window_mask)) %>%
  # Third round of masking - If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", dplyr::lag(window_mask)) | grepl("Masked", dplyr::lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P4",window_mask)) %>%
  # Fourth round of masking - If one region flanking is masked and variant count is defined as an outlier, then mask the region
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", dplyr::lag(window_mask)) | grepl("Masked", dplyr::lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P5",window_mask)) %>%
  # Exentend masked regions by 10 windows
  # First go in leading direction
  dplyr::mutate(!!!leads(window_mask)) %>%
  dplyr::ungroup() %>%
  # Make long
  tidyr::gather(lead_step, mask, -(CHROM:window_mask)) %>%
  # Convert column names from leads function to numeric number of windows away
  dplyr::mutate(lead_step = as.numeric(factor(gsub("lag_window","", lead_step)))) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  dplyr::group_by(CHROM,START_BIN,END_BIN) %>%
  dplyr::distinct(mask,.keep_all=T) %>%
  # If there is only one mask flag per group, keep it. If more than one - mark as leading mask (deals with duplicates that result from gather function)
  dplyr::mutate(lead_mask = ifelse(n()==1, mask,
                                   ifelse(n() > 1, "Masked_Lead"))) %>%
  # Take distinct chrom, start bin, while maintaining the distance to the original masked region
  dplyr::arrange(CHROM, START_BIN, desc(lead_step)) %>%
  dplyr::distinct(lead_mask, .keep_all = T) %>%
  dplyr::rowwise() %>%
  # make leading mask column, either original mask or "Lead_Mask"
  dplyr::mutate(lead_mask = ifelse(grepl("Mask", window_mask), window_mask, lead_mask)) %>%
  dplyr::select(-mask) %>%
  dplyr::group_by(CHROM) %>%
  # Repeat above operation for other side of window
  dplyr::mutate(., !!!lags(window_mask)) %>%
  dplyr::ungroup() %>%
  tidyr::gather(lag_step, mask, -(CHROM:lead_mask)) %>%
  dplyr::mutate(lag_step = as.numeric(factor(gsub("lag_window","", lag_step)))) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  dplyr::group_by(CHROM,START_BIN,END_BIN) %>%
  dplyr::distinct(mask,.keep_all=T) %>%
  dplyr::mutate(lag_mask = ifelse(n()==1, mask,
                                  ifelse(n() > 1, "Masked_Lag")))%>%
  dplyr::arrange(CHROM, START_BIN, desc(lag_step)) %>%
  dplyr::distinct(lag_mask, .keep_all = T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(lag_mask = ifelse(grepl("Mask", window_mask), window_mask, lag_mask)) %>%
  dplyr::select(-mask) %>%
  dplyr::mutate(final_mask = ifelse(grepl("Mask",window_mask), window_mask,
                                    ifelse(grepl("Mask",lead_mask), lead_mask,
                                           ifelse(grepl("Mask", lag_mask), lag_mask, "Pass")))) 


# Counts
m_windows_pr <- m_windows %>%
  dplyr::mutate(final_mask_loose = ifelse((final_mask %in% c("Masked_Lead", "Masked_Lag") & COUNT < as.numeric(args[6])) | (final_mask %in% c("Masked_Outlier") & COUNT < as.numeric(args[6]))  , "Pass", final_mask))  %>%
  # dplyr::rowwise() %>%
  dplyr::group_by(CHROM,GENOMIC_REGION)%>%
  dplyr::mutate(final_mask_loose_fill = ifelse(grepl("Mask", dplyr::lag(final_mask_loose)) & grepl("Mask",dplyr::lead(final_mask_loose)) & !grepl("Mask",final_mask_loose), "Masked_Two_Flank_Outlier_2", final_mask_loose)) 

ggplot(m_windows_pr)+
  facet_grid(CHROM~., space = "free", scales = "free")+
  geom_point(aes(x = MID_BIN/1e6, y = COUNT, color = final_mask_loose_fill)) +
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -5, yend = -5, color = final_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", final_mask_loose)))+
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -10, yend = -10, color = final_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", final_mask)))+
  theme_bw(18) +
  scale_color_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))

ggplot(m_windows_pr)+
  facet_grid(CHROM~., space = "free", scales = "free")+
  geom_point(aes(x = MID_BIN/1e6, y = COUNT, color = final_mask_loose_fill),
             data = m_windows_pr %>% dplyr::filter(!grepl("Mask", final_mask_loose_fill))) +
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -5, yend = -5, color = final_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", final_mask_loose)))+
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = -10, yend = -10, color = final_mask), size = 3,
               data = m_windows_pr %>% dplyr::filter(grepl("Mask", final_mask)))+
  theme_bw(18) +
  scale_color_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))


ggplot(m_windows_pr)+
  facet_grid(final_mask_loose_fill~., scales = "free")+
  geom_histogram(aes(x = COUNT), alpha = 0.5) +
  theme_bw(18) +
  scale_fill_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Variants/kb")+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))

ggplot(m_windows_pr)+
  facet_grid(final_mask_loose_fill~., scales = "free")+
  geom_histogram(aes(x = log(COVERAGE)), alpha = 0.5) +
  theme_bw(18) +
  scale_fill_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "log(Coverage)")+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))

ggplot(m_windows_pr)+
  facet_grid(final_mask_loose_fill~., scales = "free")+
  geom_histogram(aes(x = total_SV_CT), alpha = 0.5) +
  theme_bw(18) +
  scale_fill_viridis_d(direction = -1)+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "log(Coverage)")+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray80"))

