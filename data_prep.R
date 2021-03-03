library(tidyverse)
library(RHRV)
library(lubridate)
source("utils.R")

# Reading data

files <- list.files(path = "data/24EI/cleaned", pattern = "*.txt", full.names = T)
dat <- purrr::map(files, .f = read_file)
purrr::map2(files, dat, .f = write_file)
files_1 <- list.files(path = "data/24EI/cleaned", pattern = "*-1.*_reduced_RR_interval.txt.trimmed", full.names = T)
files_5 <- list.files(path = "data/24EI/cleaned", pattern = "*-5.*_reduced_RR_interval.txt.trimmed", full.names = T)
files <- c(files_1, files_5)
names <- stringr::str_match(files, "data\\/24EI\\/cleaned\\/(.*)_reduced_RR_interval.txt.trimmed") 
dat <- purrr::map(files, .f = read_interval)
names <- as_tibble(names) %>% 
    separate(V2, into = c("sample", NA), remove = FALSE, sep = "\\.")
names(names) <- c("path", "file", "sample")
names(dat) <- names$file

# Filtering 

filter_dat <- purrr::map(dat, function(x){
    FilterNIHR(x, long = 100, last = 10, minbpm = 40, maxbpm = 200)
})

purrr::map2(filter_dat, names(filter_dat), function(x, y){
    pdf(file = paste0("data/24EI/plots/",y,".pdf"))
    PlotNIHR(x, main = paste0("Non-interpolated instantaneous heart rate: ", y))
    plot(x$Beat$RR, type = "l",
         main = paste0("RR serie: ", y),
         ylab = "RR (ms)")
    dev.off()
})

save(dat, file = "data/24EI/Rda/raw_dat.Rda")
save(filter_dat, file = "data/24EI/Rda/filter_dat.Rda")
# Monitoring number of interval removed in each edf

RR_filtered <- purrr::map2(dat, filter_dat, function(x, y) c(nrow(x$Beat), nrow(y$Beat)))
RR_filtered <- data.frame(do.call(rbind, RR_filtered))
names(RR_filtered) <- c("before", "after")
RR_filtered$removed <- RR_filtered$before - RR_filtered$after
RR_filtered$percentage <- 100*RR_filtered$removed/RR_filtered$before


# Making segments and extract features (save to Rda files)

# Need to change "-5" here
load("data/24EI/Rda/filter_dat.Rda")
tmp <- purrr::map2(filter_dat, names(filter_dat), extract_feature)

extract_hrv <- function(dat, name){
    segments <- make_segments(dat, 300, 45) # 600 90 1200 180 1800 270
    segments <- purrr::map(segments, .f = make_features, TRUE, TRUE, TRUE)
    name <- paste0(name, "-5") # -10 -20 -30
    save(segments, file = paste0("Rda", name, ".Rda"))
    return(NULL)
}


# Making feature table. (Reading Rda files corresponding to samples)
# Need to change "-5" here

rdas <- list.files(path = "data/24EI/Rda", pattern = "*-5.Rda", full.names = T)
rdas <- purrr::map(rdas, .f = function(rda){
    load(rda)
    name <- stringr::str_match(rda, "data\\/24EI\\/Rda\\/(.*)-5.Rda")[2] 
    tmp <- purrr::map(segments, hrv2tbl)
    tmp <- dplyr::bind_rows(tmp)
    tmp$id <- paste0(name,"-",seq_along(1:nrow(tmp)))
    tmp$file <- name
    tmp$sample <- ifelse(grepl("\\.", name),
                        stringr::str_sub(name, 1, 
                                  stringr::str_locate(name, "\\.")[1,1] -1), 
                        name)
    return(tmp)
})
dat_tbl <- dplyr::bind_rows(rdas)
# only for 5 mins, remove ULF_mean
dat_tbl <- dat_tbl %>% filter(!is.na(LF_mean)) %>% 
    select(-c("ULF_mean", "size"))
save(dat_tbl, file = "data/24EI/Rda/feature-5min.Rda")

# further treatment of 5-min

load("data/24EI/features/feature-5min.Rda")
load("data/24EI/clinical/cli.Rda")
class <- cli %>% 
    dplyr::select(c("sample", "usubjid", "class", "timepoint", "ablett_score")) %>% 
    unique()
dat_tbl_original <- dat_tbl
new_names <- c("adf_type1_p_value", "adf_type2_p_value", "adf_type3_p_value",
               "kpss_type1_p_value", "kpss_type2_p_value", "kpss_type3_p_value")
dat_tbl <- dat_tbl_original %>% 
    filter(RR_count >=200 & RR_count <= 1000) %>% 
    left_join(class, by = "sample") %>% 
    mutate(class = as.factor(class)) %>% 
    dplyr::select(-c("SDANN", "SDNNIDX")) %>% 
    rename_at(vars(contains("type")), ~ new_names) %>% 
    rename(psd_slope_HfGn = psd_slop_HfGn,
           psd_spectral_index_beta = psd_slope_beta) %>% 
    mutate(psd_Hurst_exponent = ifelse(is.na(psd_slope_HfGn), 
                                       psd_slope_HfBm,
                                       psd_slope_HfGn)) %>% 
  #  relocate(psd_Hurst_exponent, .after = "psd_spectral_index_beta") %>% 
    dplyr::select(-c("psd_slope_HfBm", "psd_slope_HfGn", "size")) %>% 
    filter(!is.na(ULF_mean)) 
save(dat_tbl, file = "data/24EI/features/feature-5min-model.Rda")

log_dat <- dat_tbl %>% 
    mutate(ablett_score = as_factor(ablett_score)) %>%  
    mutate_at(c(1:8, 10:18, 21:22, 25, 27:30, 34:36, 41:44, 46, 48:49), ~ log(.x + 1)) %>% 
    relocate(psd_Hurst_exponent, .after = "psd_spectral_index_beta") %>% 
    relocate(pNN50, .after = "pNN40") %>% 
    relocate(HRV_mean, .after = "HF_mean")
save(log_dat, file = "data/24EI/features/feature-5min-model-transform.Rda")

normality <- cbind(skim_log$numeric.skewness, skim_dat$numeric.skewness, 
                   skim_log$numeric.kurtosis, skim_dat$numeric.kurtosis) %>% 
    data.frame %>% 
    mutate_all(as.numeric) %>% 
    mutate(variable = skim_log$skim_variable)
colnames(normality) <- c("skewness_after", "skewness_before", "kurtosis_after", "kurtosis_before", "variable")
normality <- normality %>% 
    dplyr::mutate(skewness_change = abs(skewness_before) - abs(skewness_after),
                  kurtosis_change = abs(kurtosis_before) - abs(kurtosis_after))

# reading clinical data

cli_enr <- readxl::read_excel("data/24EI/clinical/clinical.xlsx", sheet = "ENR", col_names = TRUE)
cli_day1 <- readxl::read_excel("data/24EI/clinical/clinical.xlsx", sheet = "Day1", col_names = TRUE)
cli_day5 <- readxl::read_excel("data/24EI/clinical/clinical.xlsx", sheet = "Day5", col_names = TRUE)
cli_events <- readxl::read_excel("data/24EI/clinical/clinical.xlsx", sheet = "Events", col_names = TRUE)


tmp_1 <- cli_day1 %>% 
    separate(USUBJID, into = c(NA, "id"), remove = FALSE, sep = "-") %>% 
    mutate(class = ifelse(ABLETT_SCORE == "NOSPASM" | ABLETT_SCORE == "SPASM", 0, 1),
           sample = paste0(id, "-1"),
           timepoint = 1)
tmp_5 <- cli_day5 %>% 
    separate(USUBJID, into = c(NA, "id"), remove = FALSE, sep = "-") %>% 
    mutate(class = ifelse(ABLETT_SCORE == "NOSPASM" | ABLETT_SCORE == "SPASM", 0, 1),
           sample = paste0(id, "-5"),
           timepoint = 5)

cli_day <- full_join(rbind(tmp_1, tmp_5), names, by = c("sample")) %>% 
    select(-c("path", "id", "EVENT", "INITIALS", "SUBJID", "enteredtime", "entry",
              "STUDYID", "SITEID"))
cli <- full_join(cli_day, cli_enr, by = "USUBJID") %>% 
    select(-c("EVENT", "INITIALS", "SUBJID", "enteredtime", "entry", 
              "STUDYID", "SITEID", "HOSP_NUMBER"))
names(cli) <- stringr::str_to_lower(names(cli))

cli <-  cli %>%  separate(time, into = c(NA, "time"), sep = " ", remove = T) %>% 
    mutate(date_time = as_datetime(paste(as.character(day), 
                                         as.character(time), sep = " ")),
           ablett_score = as.numeric(factor(ablett_score, 
                                            levels = c("NOSPASM", "SPASM", "SPASMFAILURE", "SPASMANSD"))),
           timepoint = as_factor(timepoint)) %>% 
    mutate_at(vars("tracheostomy", "ventilation",
                   "ansd_hr", "ansd_hyper_sbp", "ansd_hypo_sbp", "ansd_temp","ansd_bp", "ansd",
                   "vap_temp", "vap_wbc", "vap_ps", "vap_xray", "vap", "vap_bal",
                   "uti", "uti_temp", "uti_wbc", "uti_nitrit", "uti_culture",
                   "other_infection"), 
              .fun = recode, "INO_NO" = 0, "INO_YES" = 1)

save(cli, file = "data/24EI/Rda/cli.Rda")

# clinical data sanity check
# rule 1

cli_day1$timepoint = 1
cli_day5$timepoint = 5
cli_day <- rbind(cli_day1, cli_day5)
cli_day <- cli_day %>%  separate(time, into = c(NA, "time"), sep = " ", remove = T) %>% 
                    mutate(date_time = as_datetime(paste(as.character(day), 
                                                         as.character(time), sep = " ")))  %>% 
                    mutate(ablett_score = as.numeric(factor(ablett_score, 
                                                    levels = c("NOSPASM", "SPASM", "SPASMFAILURE", "SPASMANSD")))) %>% 
                    mutate_at(vars("tracheostomy", "ventilation",
                           "ansd_hr", "ansd_hyper_sbp", "ansd_hypo_sbp", "ansd_temp","ansd_bp", "ansd",
                           "vap_temp", "vap_wbc", "vap_ps", "vap_xray", "vap", "vap_bal",
                           "uti", "uti_temp", "uti_wbc", "uti_nitrit", "uti_culture",
                           "other_infection"), 
                      .fun = recode, "INO_NO" = 0, "INO_YES" = 1) %>% 
                    mutate(sample = paste0(str_sub(usubjid, 5), "-", timepoint)) %>% 
                    select(c("sample", "date_time", 
                            "ablett_score", 
                            "tracheostomy", "trach_day",
                             "ventilation", "ventilation_dtc",
                             "ansd_dtc", "ansd_hr", "ansd_hyper_sbp", "ansd_hypo_sbp", "ansd_temp","ansd_bp", "ansd", 
                             "vap_dtc", "vap_temp", "vap_wbc", "vap_ps", "vap_xray", "vap", "vap_bac_result", "vap_bal",
                              "uti", "uti_dtc", "uti_temp", "uti_wbc", "uti_nitrit", "uti_culture", "uti_bac_result",
                                "other_infection", "infection_dtc", "other_comp"))
    
check_ablett4 <- cli_day %>% 
                    select("sample", "ablett_score", "ansd") %>% 
                    filter((ablett_score != 4 & ansd == 1) | (ablett_score == 4 & ansd != 1))

check_ablett12 <- cli_day %>% 
                    select("sample", "ablett_score", "ansd", "ventilation", "tracheostomy") %>%
                    mutate(score_sum = rowSums(across(3:5))) %>% 
                    filter((ablett_score <=2 & score_sum != 0) | (ablett_score > 2 & score_sum == 0))

check_ablett3 <- cli_day %>% 
                    select("sample", "ablett_score", "ventilation", "tracheostomy") %>% 
                    mutate(score_sum = rowSums(across(3:4))) %>%
                    filter((ablett_score == 3 & ventilation == 0 & tracheostomy == 0) |
                            ablett_score < 3 & score_sum !=0)
check_ansd <- cli_day %>% 
                select("sample", starts_with("ansd"), -c("ansd_dtc")) %>% 
                mutate(ansd_sum = rowSums(across(2:6), na.rm = T)) %>% 
                filter((ansd_sum < 3 & ansd == 1) | (ansd_sum >= 3 & ansd == 0))

check_vap1 <-  cli_day %>% 
                    select("sample", "ventilation", "tracheostomy", "vap_dtc") %>% 
                    filter(!is.na(vap_dtc) & ventilation == 0 & tracheostomy == 0)

check_vap2 <-  cli_day %>% 
                    select("sample", "ventilation", "tracheostomy", 
                           "vap_temp", "vap_wbc", "vap_ps", "vap_xray", "vap") %>% 
                    mutate(vap_sum = rowSums(across(4:7), na.rm = T)) %>% 
                    filter(((ventilation == 1 | tracheostomy == 1) & vap_sum >= 2 & vap == 0) |
                            (((ventilation == 0 & tracheostomy == 0) | vap_sum < 2) & vap == 1))
# uti more than 2
check_uti <- cli_day %>% 
                select("sample", starts_with("uti")) %>% 
                mutate(uti_sum = rowSums(across(4:7), na.rm = T)) %>% 
                filter((uti_sum < 2 & uti == 1) | (uti_sum >= 2 & uti == 0 ))

check_ablett_transition <- as.data.frame(cbind(cli_day1$USUBJID, cli_day1$ABLETT_SCORE, 
                                               cli_day5$USUBJID, cli_day5$ABLETT_SCORE)) %>% 
                            mutate(change = V2 == V4) %>% 
                            mutate(day_1 = as.numeric(factor(V2, 
                                            levels = c("NOSPASM", "SPASM", "SPASMFAILURE", "SPASMANSD"))),
                                   day_5 = as.numeric(factor(V4, 
                                            levels = c("NOSPASM", "SPASM", "SPASMFAILURE", "SPASMANSD"))),
                                   day_1_class = ifelse(day_1 <= 2, 0, 1),
                                   day_5_class = ifelse(day_5 <= 2, 0, 1),
                                   transition_ablett = day_5 - day_1,
                                   transition_class = day_5_class - day_1_class,
                                   sample = stringr::str_sub(V1,5))


## Preparation of clinical data

load(file = "data/24EI/clinical/cli.Rda")

## Prepration of drug data



# images for presentation

load("data/24EI/Rda/050-5.Rda")
seg <- segments[[1]]    
PlotPowerBand(seg, indexFreqAnalysis = 1, ymax = 800, ymaxratio = 50,
              hr = TRUE, normalized = T, 
              Indexes = "ULF")
plot_power_band(seg, indexFreqAnalysis = 1, hr = TRUE, normalized = T)
spectrogram <- PlotSpectrogram(HRVData = seg, 
                               size = 10, shift = 2,
                               scale = "logarithmic", 
                               freqRange = c(0.03, 0.4),
                               main = "Spectrogram of 050-5-1")



library(RHRV)
seg <- CreateNonLinearAnalysis(seg)
SurrogateTest(seg, oneSided = F,
              significance = 0.05, K = 5, useFunction = asymmetryStatistic, doPlot = T)








