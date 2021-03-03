library(RHRV)
library(tidyverse)
library(strucchange)
my_skim <- skimr::skim_with(numeric = skimr::sfl(skewness = ~ e1071::skewness(., na.rm = TRUE), 
                                                 kurtosis = ~ e1071::kurtosis(., na.rm = TRUE)))

pNNX <- function(RR, gaps){
  RRdiff <- abs(diff(RR))
  unlist(map(gaps, function(x, ...) 100*sum(RRdiff > x)/length(RRdiff), RRdiff))
}
read_file <- function(file){
    dat <- read.csv(file, header = F)
    dat <- dat[-1,] # remove the first interval, usually uncompleted one
    return(dat)
}
write_file <- function(file, dat){
    readr::write_delim(as.data.frame(dat), path = paste0(file,".trimmed"),
                       delim = "", quote_escape = "none", col_names = F)
}

read_interval <- function(file){
    hrv.data <-  CreateHRVData()
    hrv.data <-  SetVerbose(hrv.data, Verbose = TRUE)
    hrv.data <- LoadBeatRR(hrv.data, file, 
                           RecordPath = "/Users/haihb/Documents/Work/Oucru/innovation", 
                           scale = 1)
    hrv.data <- BuildNIHR(hrv.data)
}

make_segments <- function(dat, segment = 300, overlap = 45){
    breaks <- seq(1, max(dat$Beat$Time)/(segment-overlap))
    starts <- (breaks-1)*(segment - overlap)
    ends <- starts + segment
    segments <- purrr::map2(starts, ends, function (s, e, ...){
        segment <- ExtractTimeSegment(dat, s, e)
        # beats <- dat$Beat[which(dat$Beat$Time > s & dat$Beat$Time <= e),]
        # segment <- CreateHRVData(Verbose = TRUE)
        # segment$Beat <- beats
        return(segment)
    },
    dat)
    return(segments)
}

# Feature extraction
make_features <- function(seg, time = TRUE, frequency = TRUE, nonlinear = TRUE){
    
    if(time){
        seg <- CreateTimeAnalysis(seg)
       # remove SDANN and SDNNIDX for 5 mins
        # unlist(seg$TimeAnalysis, recursive = TRUE)
    }
    if(frequency){
        # message(paste0("treating frequency segid: ", index))
        seg <- CreateFreqAnalysis(seg)
        seg <- InterpolateNIHR(seg)
        try({
            seg <- CalculatePowerBand(seg, indexFreqAnalysis= 1,
                                      type = "wavelet", wavelet = "d4", 
                                      bandtolerance = 0.01,
                                      relative = FALSE,
                                      ULFmin = 0, ULFmax = 0.03,
                                      VLFmin = 0.03, VLFmax = 0.05, 
                                      LFmin = 0.05, LFmax = 0.15, 
                                      HFmin = 0.15, HFmax = 0.4)
            seg$FreqAnalysis[[1]]$ULF_mean <- sum(seg$FreqAnalysis[[1]]$ULF)/length(seg$FreqAnalysis[[1]]$ULF)
            seg$FreqAnalysis[[1]]$VLF_mean <- sum(seg$FreqAnalysis[[1]]$VLF)/length(seg$FreqAnalysis[[1]]$VLF)
            seg$FreqAnalysis[[1]]$LF_mean <- sum(seg$FreqAnalysis[[1]]$LF)/length(seg$FreqAnalysis[[1]]$LF)
            seg$FreqAnalysis[[1]]$HF_mean <- sum(seg$FreqAnalysis[[1]]$HF)/length(seg$FreqAnalysis[[1]]$HF)
            seg$FreqAnalysis[[1]]$LFHF_mean <- sum(seg$FreqAnalysis[[1]]$LFHF)/length(seg$FreqAnalysis[[1]]$LFHF)
            seg$FreqAnalysis[[1]]$HRV_mean <- sum(seg$FreqAnalysis[[1]]$HRV)/length(seg$FreqAnalysis[[1]]$HRV)
            seg$FreqAnalysis[[1]]$LFnu_mean <- 100*sum(seg$FreqAnalysis[[1]]$LF/
                                                           (seg$FreqAnalysis[[1]]$LF + seg$FreqAnalysis[[1]]$HF))/length(seg$FreqAnalysis[[1]]$LF)
            seg$FreqAnalysis[[1]]$HFnu_mean <- 100*sum(seg$FreqAnalysis[[1]]$HF/
                                                           (seg$FreqAnalysis[[1]]$LF + seg$FreqAnalysis[[1]]$HF))/length(seg$FreqAnalysis[[1]]$HF)
        })
        # tiff(filename = paste0("data/24EI/spectrograms", name,"_", index,".tiff"),
        #      res = 300, width = 4, height = 4, units = 'in')
        # spectrogram <- PlotSpectrogram(HRVData = seg, 
        #                                size = 30, shift = 5,
        #                                scale = "logarithmic", 
        #                                freqRange = c(0, 0.4),
        #                                showEpLegend = F,
        #                                main = "",
        #                                ylab = "",
        #                                xlab = "", axes = F)
        # dev.off()
    }
    
    if(nonlinear){
      # message(paste0("treating frequency segid: ", index))
      try({
        seg <- CreateNonLinearAnalysis(seg)
        seg <- SurrogateTest(seg,
                             indexNonLinearAnalysis = 1,
                             oneSided = T,
                             significance = 0.05, 
                             K = 1000, 
                             useFunction = asymmetryStatistic, doPlot = F)
        
        seg <- NonlinearityTests(seg, 
                          indexNonLinearAnalysis = 1)
        # Phase space reconstruction
        time.lag <- CalculateTimeLag(seg, technique = "ami",
                                     lagMax = 500, doPlot = F, units = "Bits",
                                     method = "first.e.decay")
        embedding.dim <- CalculateEmbeddingDim(seg, numberPoints = 5000,
                                               timeLag = time.lag,
                                               threshold = 0.95,
                                               maxEmbeddingDim = 20,
                                               doPlot = F)
        max.embedding.dim = embedding.dim + 5
        # takens <- BuildTakens(seg, embedding.dim, time.lag)
        # plot3D::scatter3D(takens[,1], takens[,2], takens[,3], phi = 20, theta = 30, type = "l", col = 1)
        # correlation dimention - D2 and information dimension - D1
        seg <- NonLinearNoiseReduction(seg,
                                       embeddingDim = embedding.dim,
                                       ECGsamplingFreq = 256)
        
        seg <- CalculateCorrDim(seg, timeLag = time.lag,
                                minEmbeddingDim = embedding.dim,
                                maxEmbeddingDim = max.embedding.dim,
                                minRadius = 0.001,
                                maxRadius = 100,
                                pointsRadius = 100,
                                doPlot = F,
                                corrOrder = 2)
        corr.matrix <- seg$NonLinearAnalysis[[1]]$correlation$computations$corr.matrix
        radius <- seg$NonLinearAnalysis[[1]]$correlation$computations$radius
        range <- get_range(corr.matrix, radius)
        if(length(range) == 0)
          range <- NULL
        seg <- EstimateCorrDim(seg, regressionRange = range,
                               useEmbeddings = embedding.dim:max.embedding.dim,
                               doPlot = F)
        
        
        seg <- CalculateSampleEntropy(seg,
                                     doPlot = F)
        entropy.matrix <- seg$NonLinearAnalysis[[1]]$sampleEntropy$computations$sample.entropy
        radius <- seg$NonLinearAnalysis[[1]]$sampleEntropy$computations$radius
        range <- get_range(entropy.matrix, radius)
        if(length(range) == 0)
          range <- NULL
        seg <- EstimateSampleEntropy(seg, regressionRange = range, doPlot = F)
       
        
        # seg <- CalculateInfDim(seg, timeLag = time.lag,
        #                        minEmbeddingDim = embedding.dim,
        #                        maxEmbeddingDim = max.embedding.dim,
        #                        doPlot = F,
        #                        radius = 0.01, 
        #                        numberPoints = 50)
        # infdim.matrix <- seg$NonLinearAnalysis[[1]]$infDim$computations$log.radius
        # radius <- seg$NonLinearAnalysis[[1]]$infDim$computations$fixed.mass
        # range <- get_range(infdim.matrix, radius)
        # seg <- EstimateInfDim(seg, regressionRange = c(0.004, 0.005))
        # 
        
        seg <- CalculateMaxLyapunov(seg, 
                                    minEmbeddingDim = embedding.dim,
                                    maxEmbeddingDim = max.embedding.dim,
                                    timeLag = time.lag,
                                    radius = 60,
                                    doPlot = F, 
                                    minRefPoints = 1000,
                                    numberTimeSteps = 25
                                    )
        lya.matrix <- seg$NonLinearAnalysis[[1]]$lyapunov$computations$s.function
        radius <- seg$NonLinearAnalysis[[1]]$lyapunov$computations$time
        range <- get_range(lya.matrix, radius)
        if(length(range) == 0)
          range <- NULL
        seg <- EstimateMaxLyapunov(seg, regressionRange = sort(range),
                                   useEmbeddings = embedding.dim:max.embedding.dim,
                                   doPlot = F)
        
        radius_range <- seq(1,30, 1)
        rec_percentage <- unlist(map(radius_range, function(r){
          seg <- RQA(seg, 
                     indexNonLinearAnalysis = 1,
                     embeddingDim = embedding.dim,
                     timeLag = time.lag,
                     lmin = 2, vmin = 2,
                     radius = r,
                     doPlot = F)
          rec <- seg$NonLinearAnalysis[[1]]$rqa$REC
        }))
        r <- radius_range[which(sort(c(rec_percentage, 0.2)) == 0.2) - 1] # choose r nearest to 0.2%
        seg <- RQA(seg, 
                   indexNonLinearAnalysis = 1,
                   embeddingDim = embedding.dim,
                   timeLag = time.lag,
                   lmin = 2, vmin = 2,
                   radius = r,
                   doPlot = F)
        seg <- PoincarePlot(seg, 
                            confidence = 0.95,
                            timeLag = time.lag, 
                            doPlot = F)
        seg <- CalculateDFA(seg,
                            windowSizeRange = c(5, 200),
                            npoints = 100,
                            doPlot = T)
        df <- data.frame("f" = log(seg$NonLinearAnalysis[[1]]$dfa$computations$fluctuation.function),
                         "t" = log(seg$NonLinearAnalysis[[1]]$dfa$computations$window.sizes))
        # find number of breakpoints on log-log plot. breakpoints in terms of 
        # slope and intercept - minimum BIC
        fit_cp <- strucchange::breakpoints(f ~ 1 + t, data = df, h = 10)
        cp <- c(1, fit_cp$breakpoints, length(df$t))
        
        map(c(2:length(cp)), function(i,...){
          range <- c(df$t[cp[i-1]], df$t[cp[i]])
          range <- exp(1)^range
          seg <<- EstimateDFA(seg, regressionRange = range)
        }, seg, df, cp)
        
        seg <- CalculatePSD(seg,
                            indexFreqAnalysis = 1,
                            method = "pgram",
                            doPlot = F)
        
        seg <- EstimatePSDSlope(seg, 
                                indexFreqAnalysis = 1,
                                indexNonLinearAnalysis = 1,
                                regressionRange = c(0.005, 0.03), # estimated at low freq
                                doPlot = F)
      })
    }
    seg$niHR <- mean(seg$Beat$niHR)
    seg$iHR <- mean(seg$HR)
    tmp <- pNNX(seg$Beat$RR, c(10, 20, 30, 40))
    seg$TimeAnalysis[[1]]$pNN10 <- tmp[1]
    seg$TimeAnalysis[[1]]$pNN20 <- tmp[2]
    seg$TimeAnalysis[[1]]$pNN30 <- tmp[3]
    seg$TimeAnalysis[[1]]$pNN40 <- tmp[4]
    seg$StationaryTest <- stationary_test(seg)
    
    return(seg)
}

asymmetryStatistic <- function(x){
    x.len = length(x)
    mean(x[1:(x.len - 1)] * x[2:(x.len)]^2 - x[1:(x.len - 1)]^2 * x[2:(x.len)])
}

pNNX <- function(RR, gaps){
    RRdiff <- abs(diff(RR))
    unlist(map(gaps, function(x, ...) 100*sum(RRdiff > x)/length(RRdiff), RRdiff))
}
# Unit root test for stationarity
# Null hypothesis of adf ...
# Null hypothesis of kpss ...

stationary_test <- function(seg){
      adf <- aTSA::adf.test(seg$Beat$RR)
      kpss <- aTSA::kpss.test(seg$Beat$RR)
      adf.df <- data.frame("adf p-value" = c(adf$type1[1,3], adf$type2[1,3], adf$type3[1,3]))
      kpss.df <- data.frame("kpss p-value" = kpss[,3])
      stationary_stat <- c(t(adf.df), t(kpss.df))
      names(stationary_stat) <- c("adf_type1", "adf_type2", "adf_type3", 
                                  "kpss_type1", "kpss_type2", "kpss_type3")
    # stationary_stat <- as.data.frame(matrix(data = unlist(stationary.data, recursive = T), ncol = 6, byrow = T))
    # colnames(stationary_stat) <- c("adf1", "adf2", "adf3", "kpss1", "kpss2", "kpss3")
    
    # stationary_stat <- as_tibble(stationary_stat)
    # x <- stationary_stat %>% 
    #       mutate("adf1_call" = ifelse(adf1 >=0.05, "non-stationary", "stationary"),
    #       "kpss1_call" = ifelse(kpss1 >=0.05, "stationary with deterministic trend", "non-stationary"),
    #       "adf2_call" = ifelse(adf2 >=0.05, "non-stationary", "stationary"),
    #       "kpss2_call" = ifelse(kpss2 >=0.05, "stationary with deterministic trend", "non-stationary"),
    #       "adf3_call" = ifelse(adf3 >=0.05, "non-stationary", "stationary"),
    #       "kpss3_call" = ifelse(kpss3 >=0.05, "stationary with deterministic trend", "non-stationary"))
    return(as_tibble_row(stationary_stat))
}

# Null hypothesis: Data comes from a linear stochastic process
nonlinear_test <- function(segments, pvalue){
    nonlinear.data <- map(segments, function(seg, pvalue){
        seg <- RHRV::CreateNonLinearAnalysis(seg)
        seg <- RHRV::SurrogateTest(seg, oneSided = F,
                                  significance = pvalue, K = 100, 
                                  useFunction = asymmetryStatistic, 
                                  doPlot = F)
        test <- abs(quantile(seg$NonLinearAnalysis[[1]]$SurrogateTest$surrogates.statistics,
                 c(pvalue, 1 - pvalue))) < seg$NonLinearAnalysis[[1]]$SurrogateTest$data.statistic
        if(sum(test) > 0) 
            linear <- "rejected"
        else 
            linear <- "accepted"
        return(linear)
    }, pvalue)
    return(unlist(nonlinear.data))
}


## get_region fits lm to segments in the plot and calculate the diff in slopes and intercept of 
## lm regression across embedding dimension. These diffs are used (in get_range) to choose the segment with
## most parallel regresions across dimensions and favour the segment where intercepts are separated.
## check RHRV book for the reason.

get_region <- function(corr.matrix, radius, segnum, ancova = F, lm = T){
    names(radius) <- seq(1:length(radius))
    dat <- cbind("dim" = rownames(corr.matrix), as.data.frame(corr.matrix)) %>% 
        pivot_longer(cols = 2:ncol(.), names_to = "x", values_to = "y") %>% 
        mutate(x = as.numeric(x), dim = as_factor(dim))
    
    # pairs <- map(as.numeric(names(radius)[-c((length(radius) - segnum + 1): length(radius))]), 
    #              function(x) seq_along(x:(x + segnum)) + x - 1)
    
    pairnum <- length(radius) - segnum + 1
    pairs <- map(seq_along(1:pairnum), function(i){
      c(i: (i + segnum - 1))
    })
    pair_results <- purrr::map(pairs, function(p, ...){
        r <- radius[names(radius) %in% p]
        dat <- dat %>% filter(round(x,3) %in% round(r,3)) %>% 
            mutate(x = scales::rescale(x,  to = c(0,1)), y = y - min(y))
        if(ancova){
            fit <- glm(y ~ x*dim, data = dat)
            tmp <- anova(fit, test = "LRT")   
           }
        if(lm){
            tmp <- dat %>% 
                group_by(dim) %>% 
                do(broom::tidy(lm(y ~ x, data = .))) %>% 
                pivot_wider(id_cols = c("dim", "term"), 
                            names_from = "term",
                            values_from = c("estimate", "p.value")) %>% 
                ungroup() %>% 
                filter(across(starts_with("p.value_x")) <= 0.05) # take account into only slope
            
             tmp <- data.frame("mean_diff_slope" = abs(mean(diff(tmp$estimate_x))),
                               "mean_diff_intercept" = abs(mean(diff(tmp$'estimate_(Intercept)'))))
        }
        return(tmp)
    }, dat, radius, ancova, lm)
    names(pair_results) <- seq_along(pairs)
    regions <- map(pairs, function(pair, ...){
        data.frame("start" = radius[min(pair)], "end" = radius[max(pair)])
    }, radius)
    
    tmp <- dplyr::bind_rows(pair_results, .id = "column_label")
    
    if(ancova){
        tmp <- bind_cols("type" = rownames(tmp), tmp) %>% 
                mutate(type = stringr::str_remove(type,"\\.\\.\\.[0-9]*")) %>% 
                filter(type == "x:dim") %>% 
                pivot_wider(id_cols = c("type", "column_label"), names_from = "type", 
                        values_from = c("Resid. Dev", "Pr(>Chi)")) %>% 
            inner_join(dplyr::bind_rows(regions, .id = "column_label"), by = "column_label")
        
        idx <- as.numeric(tmp$column_label[which.max(tmp$'Pr(>Chi)_x:dim')])
        return(as_tibble(c("segnum" = segnum, 
                               "start" = tmp$start[idx],
                               "end" = tmp$end[idx], 
                               "pvalue" = max(tmp$'Pr(>Chi)_x:dim'))))
        
    }
    if(lm) {
        tmp <- tmp %>% 
            filter(!is.nan(mean_diff_slope)) %>% 
            inner_join(dplyr::bind_rows(regions, .id = "column_label"), by = "column_label") %>% 
            mutate(segnum = segnum)
        # choose idx with lowest mean diff slope and highest intercept gap
        idx <- which.max(abs(tmp$mean_diff_intercept - tmp$mean_diff_slope))
        return(tmp)
        
    }
  
}

get_range <- function(corr.matrix, radius, ancova = F, lm = T, thres = 0.1){
   regions <- map(c(3:length(radius)), # 3 points = 2 segments on plots
                   function(segnum, ...){
                       get_region(corr.matrix, radius, segnum, ancova, lm)
                   }, 
                   corr.matrix, radius)
    regions <- bind_rows(regions, .id = "column_name")
    thres <- quantile(regions$mean_diff_slope, probs = thres) # lowest 10% of slope diff
    regions <- regions %>% 
        filter(mean_diff_slope <= thres) 
    # find point that maximizes intercept diff and in 10% smallest slope diff
    id <- which.max(regions$mean_diff_intercept)
    range <- c(regions$end[id], regions$start[id])
    return(range)
}

hrv2tbl <- function(i, features, name, type){
  message(paste0("treating: ", name, " id ", i))
  seg <- features[[i]]
  
  if(length(seg$TimeAnalysis[[1]]) != 0)
      time <- unlist(seg$TimeAnalysis[[1]], recursive = T)
  else{
    time <- rep(NA, 15)
    names(time) <- c("size", "SDNN", "SDANN", "SDNNIDX", "pNN50",
                     "SDSD", "rMSSD", "IRRR", "MADRR", "TINN",
                     "HRVi", "pNN10", "pNN20", "pNN30", "pNN40")
  } 
  
  
  freq <- rep(NA, 8)

  if(length(seg$FreqAnalysis[[1]]) > 1){
    freq <- unlist(seg$FreqAnalysis[[1]][which(grepl("mean",
                                                     names(seg$FreqAnalysis[[1]])))])
  }

  names(freq) <- c("ULF_mean", "VLF_mean", "LF_mean", "HF_mean",
                   "LFHF_mean", "HRV_mean", "LFnu_mean", "HFnu_mean")

  if(length(seg$NonLinearAnalysis) > 1){
    nonlinear_out <- seg$NonLinearAnalysis[[1]]
  
  if(!is.null(nonlinear_out$correlation)) 
    correlation_dim <- c("correlation_dim_D2" = nonlinear_out$correlation$statistic)
  else
    correlation_dim <- c("correlation_dim_D2" = NA)
  
  if(!is.null(nonlinear_out$sampleEntropy)) 
      sample_entropy <- c("sample_entropy" = mean(nonlinear_out$sampleEntropy$statistic))
  else
      sample_entropy <- c("sample_entropy" = NA)
  
  if(!is.null(nonlinear_out$lyapunov)) 
    max_lyapunov_exponent <- c("max_lypunov_exponent" = nonlinear_out$lyapunov$statistic)
  else
    max_lyapunov_exponent <- c("max_lypunov_exponent" = NA)
  
  if(!is.null(nonlinear_out$rqa))
      rqa <- unlist(null_to_na(nonlinear_out$rqa[2:12]))
  else
      rqa <- rep(NA, 11)
  
  names(rqa) <- c("rqa_REC","rqa_RATIO", "rqa_DET", "rqa_DIV", "rqa_Lmax", "rqa_Lmean", 
                  "rqa_LmeanWithoutMain", "rqa_ENTR", "rqa_TREND", "rqa_LAM", "rqa_Vmax")
  if(!is.null(nonlinear_out$PoincarePlot))
      poincare <- unlist(nonlinear_out$PoincarePlot)
  else
      poincare <- c(NA, NA)
  
  names(poincare) <- c("poincare_SD1", "poincare_SD2")
  dfas <- nonlinear_out$dfa$statistic
  dfa <- c(NA, NA)
  # getting first 2 scaling exponents
  dfa[1] <- ifelse(!is.null(dfas[[1]]), dfas[[1]]$estimate, NA)
  dfa[2] <- ifelse(!is.null(dfas[[2]]), dfas[[2]]$estimate, NA)
  names(dfa) <- c("dfa_alpha1", "dfa_alpha2")
  
  if(!is.null(nonlinear_out$PSDSlope))
      psd <- unlist(null_to_na(nonlinear_out$PSDSlope$statistic))
  else
      psd <- rep(NA, 3)
  
  names(psd) <- c("psd_spectral_index_beta", "psd_slope_HfBm", "psd_slope_HfGn")
  
  surrogate_f <- ecdf(nonlinear_out$SurrogateTest$surrogates.statistics)
  surrogate_test <- c("surrogate_p_value" = surrogate_f(nonlinear_out$SurrogateTest$data.statistic))
  
  nonlinear_test <- unlist(nonlinear_out$NonlinearityTests, recursive = T) %>% 
    as_tibble_row %>% 
    mutate('McLeodLi.p.value' = mean(as.numeric(c_across(dplyr::starts_with("McLeodLi.p.values"))))) %>% 
    dplyr::select(ends_with("p.value")) %>% 
    dplyr::mutate_all(as.numeric) %>%
    dplyr::rename_all(str_to_lower) %>% 
    dplyr::rename_all(str_replace_all, "\\.", "_") %>% 
    unlist()
  stationary_test <- unlist(seg$StationaryTest)
  }
  hr <- c("iHR" = seg$iHR, "niHR" = seg$niHR)
  RR_count <- c("RR_count" = nrow(seg$Beat))
  
  # return(c(time, freq,
  #          correlation_dim, sample_entropy, max_lyapunov_exponent,
  #          rqa, poincare, dfa, psd,
  #          surrogate_test, nonlinear_test, stationary_test,
  #          hr, RR_count))

  return(c("sample" = name, "id" = i, "type" = type, time, freq,
           hr, RR_count))
}


## get_dfa_range fits a series of lm regressions in the log(fluctuation) vs log(time) plot
## to find scaling regions. The best fits are defined by the largest R squared and tie-break by
## the largest number of points (longest segment). 2 crossover regions are selected for short and long
## term (in time) fluctuation behaviors. Lowest number of point per segment theta is 10. 
# M (total number of points) is 100, estimated from RHRV.
# Cite (A criterion for the determination of optimal scaling ranges in DFA and MF-DFA)

get_dfa_range <- function(f, t, theta = 10){
  f <- log(f)
  t <- log(t)
  regions <- map(c(theta:length(f)), function(segnum,...){
    pairnum <- length(f) - segnum + 1
    pairs <- map(seq_along(1:pairnum), function(i){
      c(i: (i + segnum - 1))
    })
    r_squared <- map(seq_along(pairs), function(i, ...){
      pair <- pairs[[i]]
      r_squared <- summary(lm(f[pair] ~ t[pair]))$r.squared
      names(r_squared) <- as.character(i)
      return(r_squared)
    }, pairs, f, t)
    names(r_squared) <- rep(as.character(segnum), length(r_squared))
    return(r_squared)
  }, f, t)
  regions <- sort(unlist(regions, recursive = T), decreasing = T)[1:2]
  tmp <- map(names(regions), function(region,...){
    
  }, t)
}

# Testing segmented regression
# lm.br(f~t)
# my.lm <- lm(f~t)
# summary(my.lm)
# plot(t,f)
# my.seg <- segmented(my.lm, psi = median(t))
# my.seg1 <- segmented(my.lm, npsi = 3)
# my.seg2 <- segmented(my.lm, npsi = 3)
# my.seg3 <- segmented(my.lm, npsi = 2)
# my.seg4 <- segmented(my.lm, npsi = 5)
# summary(my.seg)
# my.fitted <- fitted(my.seg)
# my.model <- data.frame("t" = t, "fitted" = my.fitted)
# 
# ggplot(my.model, aes(x = t, y = fitted)) + geom_line() +
#   geom_point(aes(x = t, y = f), colour = "red") +
#   # geom_vline(xintercept = my.seg$psi[,2], linetype = "dashed", colour = "green") +
#   #   geom_vline(xintercept = my.seg1$psi[,2], linetype = "dashed", colour = "blue")  +
#   geom_vline(xintercept = t[fit_in$breakpoints], linetype = "dashed", colour = "magenta") +
#   geom_vline(xintercept = t[fit$breakpoints], linetype = "dashed", colour = "yellow") 


# Set of plots: divide by class (2 or 4 class), by day

plot_set3 <- function(dat, class){
  if(class == 2){
    print(dat %>% 
            # compare_means(value ~ class, data = ., method = "wilcox.test", paired = FALSE,
            #               group.by = "variable")
            ggplot(aes(x = variable, y = value, fill = class)) + geom_boxplot() +
            scale_fill_manual(name="Class",
                              labels=c("Mild","Severe"),
                              values= c("#e2cc92","#c74e1a")) +
            # stat_summary(fun = median, geom="point", size=1, color="black") +
            # stat_compare_means(method = "wilcox.test", paired = FALSE, 
            # #                  aes(label = ..p.signif..),
            # label.x = 1.5, label.y = 13) +
            ggtitle("Distribution of time features by class") + 
            ylab("value (log scale)") + theme(plot.title = element_text(hjust = 0.5))
    )
    print(dat %>% 
            # compare_means(value ~ class, data = ., method = "wilcox.test", paired = FALSE,
            #               group.by = "variable")
            ggplot(aes(x = variable, y = value, fill = class)) + geom_boxplot() +
            scale_fill_manual(name="Class",
                              labels=c("Mild","Severe"),
                              values= c("#e2cc92","#c74e1a")) +
            # stat_summary(fun = median, geom="point", size=1, color="black") +
            # stat_compare_means(method = "wilcox.test", paired = FALSE, 
            # #                  aes(label = ..p.signif..),
            # label.x = 1.5, label.y = 13) +
            ggtitle("Distribution of time features by class") + 
            ylab("value (log scale)") + theme(plot.title = element_text(hjust = 0.5)) +
            facet_wrap(vars(timepoint), 
                       labeller = labeller(timepoint = c("1" = "Day 1", "5" = "Day 5")))
    )
    
    print(dat %>% 
            # compare_means(value ~ class, data = ., method = "wilcox.test", paired = FALSE,
            #               group.by = "variable")
            ggplot(aes(x = variable, y = value, fill = timepoint)) + geom_boxplot() +
            scale_fill_manual(name="Timepoint",
                              labels=c("Day 1","Day 5"),
                              values= c("#97c4f5","#e5bfab")) +
            # stat_summary(fun = median, geom="point", size=1, color="black") +
            # stat_compare_means(method = "wilcox.test", paired = FALSE, 
            # #                  aes(label = ..p.signif..),
            # label.x = 1.5, label.y = 13) +
            ggtitle("Distribution of time features by class") + 
            ylab("value (log scale)") + theme(plot.title = element_text(hjust = 0.5)) +
            facet_wrap(vars(class), labeller = labeller(class = c("0" = "Mild", "1" = "Severe")))
    )
  }
  if(class == 4){
    print(dat %>% 
            # compare_means(value ~ class, data = ., method = "wilcox.test", paired = FALSE,
            #               group.by = "variable")
            ggplot(aes(x = variable, y = value, fill = as_factor(ablett_score))) + geom_boxplot() +
            scale_fill_manual(name="Ablett score",
                              labels=c("1","2", "3", "4"),
                              values=c("#e2cc92", "#d1a650","#da7a24","#c74e1a")) +
            # stat_summary(fun = median, geom="point", size=1, color="black") +
            # stat_compare_means(method = "wilcox.test", paired = FALSE, 
            #                  aes(label = ..p.signif..),
            # label.x = 1.5, label.y = 13) +
            ggtitle("Distribution of time features by Ablett score") + 
            ylab("value (log scale)") + theme(plot.title = element_text(hjust = 0.5))  
    )
    print(dat %>% 
            # compare_means(value ~ class, data = ., method = "wilcox.test", paired = FALSE,
            #               group.by = "variable")
            ggplot(aes(x = variable, y = value, fill = as_factor(ablett_score))) + geom_boxplot() +
            scale_fill_manual(name="Ablett score",
                              labels=c("1","2", "3", "4"),
                              values=c("#e2cc92", "#d1a650","#da7a24","#c74e1a")) +
            # stat_summary(fun = median, geom="point", size=1, color="black") +
            # stat_compare_means(method = "wilcox.test", paired = FALSE, 
            #                  aes(label = ..p.signif..),
            # label.x = 1.5, label.y = 13) +
            ggtitle("Distribution of time features by Ablett score") + 
            ylab("value (log scale)") + theme(plot.title = element_text(hjust = 0.5)) +
            facet_wrap(vars(timepoint), 
                       labeller = labeller(timepoint = c("1" = "Day 1", "5" = "Day 5")))
    )
    print(dat %>% 
            # compare_means(value ~ class, data = ., method = "wilcox.test", paired = FALSE,
            #               group.by = "variable")
            ggplot(aes(x = variable, y = value, fill = timepoint)) + geom_boxplot() +
            scale_fill_manual(name="Timepoint",
                              labels=c("Day 1","Day 5"),
                              values= c("#97c4f5","#e5bfab")) +
            # stat_summary(fun = median, geom="point", size=1, color="black") +
            # stat_compare_means(method = "wilcox.test", paired = FALSE, 
            # #                  aes(label = ..p.signif..),
            # label.x = 1.5, label.y = 13) +
            ggtitle("Distribution of time features by class") + 
            ylab("value (log scale)") + theme(plot.title = element_text(hjust = 0.5)) +
            facet_wrap(vars(ablett_score), labeller = labeller(ablett_score = c("1" = "Score 1", "2" = "Score 2",
                                                                                "3" = "Score 3", "4" = "Score 4")))
    )
  }
}

null_to_na <- function(l){
  map(l, function(x){
    if(is.null(x))
      x <- NA
    return(x)
  })
}

