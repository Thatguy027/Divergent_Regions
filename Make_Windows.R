library(tidyverse)
library(zoo)
library(rlang)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

leads <- function(var, n=10){
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(dplyr::lead(!!var, !!.x)) ) %>% 
    set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
  
}

lags <- function(var, n=10){
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(dplyr::lag(!!var, !!.x)) ) %>% 
    set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
  
}

outliers <- readr::read_tsv("CB4856_Processed_Outliers.tsv")

# Coverage - remove coverage outliers

m_windows <- outliers %>%
  dplyr::mutate(window_mask = ifelse(cov_direction %in% c("Down", "Up") | count_rank < 100, "Masked", "Pass")) %>%
  dplyr::mutate(window_mask = ifelse(window_mask == "Masked", "Masked", "Pass")) %>%
  dplyr::mutate(window_mask = ifelse(count_outlier == "Yes" & fraction_SV_bases > 0.5, "Masked_SV", window_mask)) %>%
  dplyr::mutate(window_mask = ifelse((lag(window_mask) == "Masked" & lead(window_mask) == "Masked") & !grepl("Masked", window_mask), "Masked_P1", window_mask)) %>%
  dplyr::mutate(window_mask = ifelse((lag(window_mask) == "Masked" | lead(window_mask) == "Masked") & count_outlier == "Yes","Masked_P2", window_mask)) %>%
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", lag(window_mask)) | grepl("Masked", lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P3",window_mask)) %>%
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", lag(window_mask)) | grepl("Masked", lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P4",window_mask)) %>%
  dplyr::mutate(window_mask = ifelse(((grepl("Masked", lag(window_mask)) | grepl("Masked", lead(window_mask))) & count_outlier == "Yes" & !grepl("Masked", window_mask)),"Masked_P5",window_mask)) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(!!!leads(window_mask, 10)) %>%
  dplyr::ungroup() %>%
  tidyr::gather(lead_step, mask, -(CHROM:window_mask)) %>%
  dplyr::mutate(lead_step = as.numeric(factor(gsub("lag_window","", lead_step)))) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  dplyr::group_by(CHROM,START_BIN,END_BIN) %>%
  dplyr::distinct(mask,.keep_all=T) %>%
  dplyr::mutate(lead_mask = ifelse(n()==1, mask,
                                   ifelse(n() > 1, "lead_Mask"))) %>%
  dplyr::arrange(CHROM, START_BIN, desc(lead_step)) %>%
  dplyr::distinct(lead_mask, .keep_all = T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(lead_mask = ifelse(grepl("Mask", window_mask), window_mask, lead_mask)) %>%
  dplyr::select(-mask) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(., !!!lags(window_mask, 10)) %>%
  dplyr::ungroup() %>%
  tidyr::gather(lag_step, mask, -(CHROM:lead_mask)) %>%
  dplyr::mutate(lag_step = as.numeric(factor(gsub("lag_window","", lag_step)))) %>%
  dplyr::arrange(CHROM, START_BIN) %>%
  dplyr::group_by(CHROM,START_BIN,END_BIN) %>%
  dplyr::distinct(mask,.keep_all=T) %>%
  dplyr::mutate(lag_mask = ifelse(n()==1, mask,
                                   ifelse(n() > 1, "lag_Mask")))%>%
  dplyr::arrange(CHROM, START_BIN, desc(lag_step)) %>%
  dplyr::distinct(lag_mask, .keep_all = T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(lag_mask = ifelse(grepl("Mask", window_mask), window_mask, lag_mask)) %>%
  dplyr::select(-mask) %>%
  dplyr::mutate(final_mask = ifelse(grepl("Mask",window_mask), window_mask,
                                    ifelse(grepl("Mask",lead_mask), lead_mask,
                                           ifelse(grepl("Mask", lag_mask), lag_mask, "Pass"))))
  
  
d <- data_frame(x = seq_len(100))
mutate( d, !!!lags(x, 3) )
  
  dplyr::filter(n()==1 | grepl("Masked", mask)) %>%
  dplyr::mutate(lead_mask = ifelse(grepl("Mask", window), window, 
                                   ifelse(grepl("Mask", mask), "lead_Masked", "Pass")))

test <- m_windows %>%
  dplyr::filter(n()==1)

unique(test$mask)


  dplyr::mutate(lead_window = ifelse(grepl("Masked", unique(mask)), "Lead_Mask", "Pass"))

m_windows %>% dplyr::group_by(final_mask) %>% dplyr::summarise(sum_ct = sum(COUNT))



# Counts
ggplot()+
  facet_grid(CHROM~., space = "free", scales = "free")+
  geom_line(aes(y=smoothed_count_outlier,x = MID_BIN/1e6), color = "red", data = m_windows) +
  geom_segment(aes(x = START_BIN/1e6, xend = END_BIN/1e6, y = 0, yend = 0, color = final_mask), size = 3,
               data = m_windows %>% dplyr::filter(grepl("Mask", final_mask)))+
  theme_bw(18) +
  scale_color_viridis_d()+
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count")+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray70")) + xlim(c(0,3))


