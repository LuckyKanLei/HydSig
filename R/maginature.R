data_Dresden <- read.csv("Region193_Dresden.csv")
data_Dresden$date <- as.Date(data_Dresden$date)
library(xts)
library(purrr)
library(plyr)
library(zoo)
library(rmatio)
library(HMtools)

data29 <- read.mat("33029_daily.mat")
data29[c(1,6)] <- NULL
df_data29 <- as.data.frame(data29)

data20 <- read.mat("39020_daily.mat")
data20[c(1,6)] <- NULL
df_data20 <- as.data.frame(data20)

df_data29$Date <- seq(as.Date("1999-1-1"), as.Date("2008-12-31"), length.out = 3653)
xts_Data29 <- xts(df_data29[,1:4], df_data29$Date)
xts_Data20 <- xts(df_data20, df_data29$Date)

xts_Q <- cbind(xts_Data29$Q, xts_Data20$Q)
xts_P <- cbind(xts_Data29$P, xts_Data20$P)
xts_T <- cbind(xts_Data29$T, xts_Data20$T)
xts_PET <- cbind(xts_Data29$PET, xts_Data20$PET)

apply.daily




#' @title Mean
#' @description  Mean of data, e.g. streamflow, precipitation, temperature and so on.
#' @param Data [xts]-2D(time, space) Time serious Data in difficult location.
#' @param t_scale [function / NULL] |apply.daily| for daily, |apply.weekly| for weekly, |apply.monthly| for monthly,
#' |apply.quarterly| for quarterly, |apply.yearly| for yearly, |NULL| for all of the time.
#' @importFrom xts apply.daily apply.weekly apply.monthly apply.quarterly apply.yearly
#' @return [xts / array]-2D(time, space) Mean of data.
#' @example sig_Mean(xts_Q)
#' sig_Mean(xts_Q， t_scale = NULL)
sig_Mean <- function(Data, t_scale = apply.monthly){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  if(is.null(t_scale)) return(colMeans(Data, na.rm = T))
  else return(t_scale(Data, colMeans, na.rm = T))
}

#' @title  Quantile / Percentile decreasing
#' @description Quantile / Percentile of data in decreasing, e.g. streamflow, precipitation, temperature and so on.
#' The Q5 is the high flow, Q95 is the low flow.
#' For the other scale (if you want) can use apply.yearly(Data, sig_Quantile, tile = c(5, 95)) or apply.weekly apply.monthly apply.quarterly apply.yearly.
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param tile [num / vector]-range(0-100) The tiles (positions).
#' @return [xts / array]-2D(tile, space) The data in the tiles.
#' @examples sig_Quantile(xts_Q)
#' sig_Quantile(xts_Q, tile = c(10, 50, 90))
sig_Quantile <- function(Data, tile = c(5, 95)){
  if(any(tile < 0 | tile > 100)) errorCondition("Please give the tile between 0 and 100.")
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  tile <- tile / 100.
  tile_N <- length(tile)
  sort_Data <- apply(Data, 2, sort, decreasing = T, na.last = T)
  dim_Data <- dim(Data)
  time_N <- colSums(!is.na(sort_Data))
  Quan_ <- array(rep(1:dim_Data[1], dim_Data[2]), dim = dim_Data) / array(rep(time_N, each = dim_Data[1], ), dim = dim_Data)
  smaller_Quan <- bigger_Quan <- array(0, dim = c(tile_N, dim_Data[2]))

  for (i in 1:tile_N) {
    s_tile <- apply(Quan_, 2, function(x) max(which(x <= tile[i])))
    smaller_Quan[i, ] <- sort_Data[cbind(s_tile, 1:dim_Data[2])]
    b_tile <- apply(Quan_, 2, function(x) min(which(x >= tile[i])))
    bigger_Quan[i, ] <- sort_Data[cbind(b_tile, 1:dim_Data[2])]
  }
  Quantile_ <- (smaller_Quan + bigger_Quan) / 2
  colnames(Quantile_) <- colnames(Data)
  rownames(Quantile_) <- paste0("tile_", tile)
  return(Quantile_)
}

#' @title Runoff ratio
#' @description Runoff ratio of the time serise data.
#' For the other scale (if you want) can use apply.yearly or apply.weekly apply.monthly apply.quarterly apply.yearly see the example.
#' please make sure the Streamflow and Preciptation have the same unit.
#' @param Streamflow [xts / array]-2D(time, space) Time serious Streamflow in difficult location.
#' @param Preciptation [xts / array]-2D(time, space) Time serious Preciptation in difficult location.
#' @return [array]-2D(1, space) The Runoff ratio. Or [xts]-2D(sum_time, space) for other scale with apply.XXXly.
#' @example sig_Runoff_ratio(xts_Q, xts_P)
#' apply.yearly(cbind(Q = xts_Q, P = xts_P), function(QP) sig_Runoff_ratio(QP[, grep("Q", colnames(QP))], QP[, grep("P", colnames(QP))]))
sig_Runoff_ratio <- function(Streamflow, Preciptation){
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(c(Streamflow, Preciptation)))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  return(colSums(Streamflow, na.rm = T) / colSums(Preciptation, na.rm = T))
}

#' @title Base flow
#' @description Base flow with digital filter method.
#' @references
#' [Arnold.1995]
#' [eq. 1-2]
#' @param Streamflow [xts]-2D(time, space) Time serious Streamflow in difficult location.
#' @param param_filter [num]-range() Filter parameter that enables the shape of the separation to be altered.
#' @return baseflow [xts]-2D(time, space) Time serious baseflow in difficult location.
#' @example base_Flow.DFM_matrix(xts_Q)
base_Flow.DFM_matrix <- function(Streamflow, param_filter = 0.925){
  dim_Q <- dim(Streamflow)
  time_N <- dim_Q[1]
  mat_filter <- matrix(param_filter^(rep(1:time_N,time_N) - rep(1:time_N, each=time_N)), time_N)
  mat_filter[upper.tri(mat_filter)] <- 0
  if(any(is.na(Streamflow))) {
    message("The Null (NA) values in Streamflow will replace by 0, if nessary please deal with the NAs. And the results in the same place will in NA seted.")
    Streamflow0 <- replace(Streamflow, is.na(Streamflow), 0)
    Q_diff <- diff(Streamflow0) * (1 + param_filter) * 0.5
    Q_direct <- apply(Q_diff, 2, function(Q_) mat_filter %*% Q_)
    return((Streamflow - Q_direct)[-1, ])
  }
  Q_diff <- rbind(0, diff(as.matrix(Streamflow)) * (1 + param_filter) * 0.5)
  Q_direct <- apply(Q_diff, 2, function(Q_) mat_filter %*% Q_)
  return((Streamflow - Q_direct)[-1, ])
}

#' @title Base flow
#' @description Base flow with digital filter method.
#' @references
#' [Arnold.1995]
#' [eq. 1-2]
#' @param Streamflow [xts]-2D(time, space) Time serious Streamflow in difficult location.
#' @param param_filter [num]-range() Filter parameter that enables the shape of the separation to be altered.
#' @return baseflow [xts]-2D(time, space) Time serious baseflow in difficult location.
#' @example base_Flow.DFM(xts_Q)
#' base_Flow.DFM(xts_Q, passes_N = 3, negativ_step = F)
base_Flow.DFM <- function(Streamflow, param_filter = 0.925, passes_N = 1, negativ_step = T){
  ## deal with NAs
  if(any(is.na(Streamflow))) interpolat_NA_time(Streamflow)
  dim_Q <- dim(Streamflow)
  time_N <- dim_Q[1]
  Q_diff <- as.matrix(diff(Streamflow) * (1 + param_filter) * 0.5)
  Q_direct <- (Q_diff)
  Q_direct[1, ] <- Streamflow[1, ] - apply(Streamflow, 2, min)
  if(negativ_step) {
    for (i in 2:(time_N)) {
      Q_direct[i, ] <- maxSVector(0, param_filter * Q_direct[i-1, ] + Q_diff[i, ])
    }

  } else {
    for (i in 2:(time_N)) {
      Q_direct[i, ] <- param_filter * Q_direct[i-1, ] + Q_diff[i, ]
    }
    Q_direct <- maxSVector(0, Q_direct)
  }
  Q_base <- Streamflow - Q_direct
  if(passes_N > 1) for (i in 2:passes_N) {Q_base <- fct_filter.DFM(as.xts(rev(as.xts(Q_base))), param_filter, negativ_step)}
  return(Q_base)
}

fct_filter.DFM <- function(Streamflow, param_filter = 0.925, negativ_step = T){
  dim_Q <- dim(Streamflow)
  time_N <- dim_Q[1]
  Q_diff <- as.matrix(diff(Streamflow) * (1 + param_filter) * 0.5)
  Q_direct <- (Q_diff)
  Q_direct[1, ] <- Streamflow[1, ] - apply(Streamflow, 2, min)
  if(negativ_step) {
    for (i in 2:(time_N)) {
      Q_direct[i, ] <- maxSVector(0, param_filter * Q_direct[i-1, ] + Q_diff[i, ])
    }

  } else {
    for (i in 2:(time_N)) {
      Q_direct[i, ] <- param_filter * Q_direct[i-1, ] + Q_diff[i, ]
    }
    Q_direct <- maxSVector(0, Q_direct)
  }
  return((Streamflow - Q_direct))
}

#' @title Baseflow index
#' @description Baseflow is from digital filter method with base_Flow.DFM calculated.
#' @param Streamflow [xts]-2D(time, space) Time serious Streamflow in difficult location.
#' @param param_filter [num]-range() Filter parameter that enables the shape of the separation to be altered.
#' @param t_scale [function / NULL] |apply.daily| for daily, |apply.weekly| for weekly, |apply.monthly| for monthly,
#' |apply.quarterly| for quarterly, |apply.yearly| for yearly, |NULL| for all of the time.
#' @return Baseflow index [xts]-2D(1, space) Time serious baseflow in difficult location.
#' @example sig_Baseflow_index.DFM(xts_Q, negativ_step = F, passes_N = 3, t_scale = NULL)
#' sig_Baseflow_index.DFM(xts_Q)
sig_Baseflow_index.DFM <- function(Streamflow, param_filter = 0.925, passes_N = 1, negativ_step = T, t_scale = apply.yearly){
  Q_base <- base_Flow.DFM(Streamflow, param_filter, passes_N, negativ_step)
  if(is.null(t_scale)) mat_rate <- colSums(Q_base, na.rm = T) / colSums(Streamflow, na.rm = T)
  else mat_rate <- t_scale(cbind(B = Q_base, Q = Streamflow), function(BQ) colSums(BQ[grep("B", colnames(BQ)), ], na.rm = T) / colSums(BQ[grep("Q", colnames(BQ)), ], na.rm = T))
  return(mat_rate)
}

#' @title Reshape the time series in water-year
#' @description Reshape the time series in water-year, set the start month of the water year,
#' the number of a water year will be same as the first year-number. And the data before the start month will in the end of the time series.
#' Please make sure the time series begin with the 1. Jan..
#' @param Data [xts]-2D(time, space) Time serious Data in difficult location.
#' @param wateryear_start_m [num]-|2|3|4|5|6|7|8|9|10|11|12| for Feb. ... Dec.
#' @return By water year reshaped time serious Data [xts]-2D(time, space)
#' @example sig_Baseflow_index.DFM(xts_Q, 11)
water_year <- function(Data, wateryear_start_m){
  start_y <- .indexyear(Data)[1]
  start_yd <- .indexyday(Data)[1]
  N_preData <- sum(.indexyear(Data) == start_y & .indexmon(Data) < (wateryear_start_m - 1))
  if(start_yd != 0) errorCondition("If you want to use the \"water year\" start with other month, please make sure that the data series begin with 1.Jan..")
  Data[,] <- rbind(as.matrix(Data[-(1:N_preData), ]), as.matrix(Data[(1:N_preData), ]))
  return(Data)
}


#' @param Streamflow [xts]-2D(time, space) Time serious Streamflow in difficult location.
#' @param Preciptation [xts]-2D(time, space) Time serious Preciptation in difficult location.
#' @param t_scale [function] |apply.quarterly| for quarterly, |apply.yearly| for yearly.
#' @param wateryear_start_m [num]-|2|3|4|5|6|7|8|9|10|11|12| for Feb. ... Dec.
#' @return stream elasticity [xts]-2D(time, space) in time scale
#' @example sig_Stream_elas.Sawicz(xts_Q, xts_P)
#' median(sig_Stream_elas.Sawicz(xts_Q, xts_P, wateryear_start_m = 10))
sig_Stream_elas.Sawicz <- function(Streamflow, Preciptation, t_scale = apply.yearly, wateryear_start_m = 1) {
  if(t_scale != apply.yearly && wateryear_start_m != 1) message("Water year is only for yearly availability.")
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Preciptation))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  if(wateryear_start_m > 1) {
    Streamflow <- water_year(Streamflow, wateryear_start_m)
    Preciptation <- water_year(Preciptation, wateryear_start_m)
    }
  Q_mean_t_scale <- t_scale(Streamflow, colMeans, na.rm = T)
  P_mean_t_scale <- t_scale(Preciptation, colMeans, na.rm = T)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  dQ <- diff(Q_mean_t_scale)
  dP <- diff(P_mean_t_scale)
  elas_ <- dQ / dP * rep(colMeans(P_mean_t_scale), each = dim_scale[1]) / rep(colMeans(Q_mean_t_scale), dim_scale[1])
  return(elas_[-1, ])
}


sig_Stream_elas.Sanka <- function(Streamflow, Preciptation, t_scale = apply.yearly, wateryear_start_m = 1) {
  if(t_scale != apply.yearly && wateryear_start_m != 1) message("Water year is only for yearly availability.")
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Preciptation))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  Q_mean_t_scale <- t_scale(Streamflow, colMeans, na.rm = T)
  P_mean_t_scale <- t_scale(Preciptation, colMeans, na.rm = T)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  Q_mean <- rep(colMeans(Q_mean_t_scale), each = dim_scale[1])
  P_mean <- rep(colMeans(P_mean_t_scale), each = dim_scale[1])
  elas_ <- (Q_mean_t_scale - Q_mean) / (P_mean_t_scale - P_mean) * P_mean / Q_mean
  return(elas_[-1,])
}

sig_Frequency_threshold <- function(Data, higher = T, threshold = c(1, 3, 9), times_mean = T){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")

  time_N <- dim(Data)[1]
  space_N <- dim(Data)[2]
  if(times_mean) {
    threshold_N <- length(threshold)
    threshold <- matrix(rep(apply(Data, 2, mean), each = threshold_N) * rep(threshold, space_N), threshold_N, space_N)
  } else {
    if(is.null(dim(threshold))) {
      message("All stations (or grids, the space points) will use the same threshold(s).
               If want to caculate the different thresholds for each stations please use array() / matrix() to set the \'threshold\'.")
      threshold_N <- length(threshold)
      threshold <- matrix(rep(threshold, space_N), threshold_N, space_N)
    } else {
      if(dim(threshold)[2] != space_N) errorCondition("The second dimesion of \'threshold\' must gleich as the \'Data\'.")
      else threshold_N <- dim(threshold)[1]}
  }

  dim_Data <- dim(Data)
  mat_freq <- matrix(0, threshold_N, space_N)
  colnames(mat_freq) <- colnames(Data)
  if(higher){
    for (i in 1:threshold_N) {
      mat_Diff <- Data - rep(threshold[i,], each = time_N)
      mat_freq[i,] <- colSums(mat_Diff > 0, na.rm = T) / colSums(!is.na(Data), na.rm = T)
    }
  } else {
    for (i in 1:threshold_N) {
      mat_Diff <- Data - rep(threshold[i,], each = time_N)
      mat_freq[i,] <- colSums(mat_Diff <= 0, na.rm = T) / colSums(!is.na(Data), na.rm = T)
    }
  }
  return(mat_freq)
}
sig_Duration_threshold <- function(Data, higher = T, threshold = c(1, 3, 9), times_mean = T, dur_type = mean){
  time_N <- dim(Data)[1]
  space_N <- dim(Data)[2]
  if(times_mean) {
    threshold_N <- length(threshold)
    threshold <- matrix(rep(apply(Data, 2, mean), each = threshold_N) * rep(threshold, space_N), threshold_N, space_N)
  } else {
    if(is.null(dim(threshold))) {
      message("All stations (or grids, the space points) will use the same threshold(s).
               If want to caculate the different thresholds for each stations please use array() / matrix() to set the \'threshold\'.")
      threshold_N <- length(threshold)
      threshold <- matrix(rep(threshold, space_N), threshold_N, space_N)
    } else {
      if(dim(threshold)[2] != space_N) errorCondition("The second dimesion of \'threshold\' must gleich as the \'Data\'.")
      else threshold_N <- dim(threshold)[1]}
  }

  dim_Data <- dim(Data)
  mat_dur <- matrix(0, threshold_N, space_N)
  colnames(mat_dur) <- colnames(Data)
  if(higher){
    for (i in 1:threshold_N) {
      mat_Diff <- Data - rep(threshold[i,], each = time_N) > 0
      mat_dur[i,] <- apply(mat_Diff, 2, function(D) {
        rle_Diff <- rle(as.vector(D))
        return(dur_type(rle_Diff$lengths[which(rle_Diff$values)]))})
    }
  } else {
    for (i in 1:threshold_N) {
      mat_Diff <- Data - rep(threshold[i,], each = time_N) <= 0
      mat_dur[i,] <- apply(mat_Diff, 2, function(D) {
        rle_Diff <- rle(as.vector(D))
        return(dur_type(rle_Diff$lengths[which(rle_Diff$values)]))})
    }
  }
  return(mat_dur)
}

sig_Snowday_ratio <- function(Preciptation, Temperatur, snow_threshpld = 2){
  dim_T <- dim(Temperatur)
  dim_P <- dim(Preciptation)
  if(dim_T[1] != dim_P[1] | dim_T[2] != dim_P[2]) errorCondition("Please make sure Temperatur and Preciptation have the same dimsion.")
  if(any(is.na(Preciptation)) | any(is.na(Temperatur))) message("The NA in Temperatur and Preciptation will be droped, if nessary please deal with the NAs.")

  judge_P <- Preciptation > 0
  judge_S <- judge_P & (Temperatur < snow_threshpld)
  S_ration <- colSums(judge_S, na.rm = T) / colSums(judge_P, na.rm = T)
  S_ration[is.na(S_ration)] <- 0
  return(S_ration)
}

sig_Autocor <- function(Data, time_lag = 1){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  time_N <- dim(Data)[1]
  return(apply(Data, 2, function(D) {
    D <- as.numeric(D[which(!is.na(D))])
    cor(D[-((length(D) - time_lag + 1):length(D))], D[-(1:time_lag)])}))
}



sig_Rising_limb_density <- function(Data, decrease_toleranz = 0, min_rising_duration = 1){
  ## deal with NAs
  if(any(is.na(Data))){
    message("The NA in \'Data\' will be linear interpolated with the last and next data, that is not NA, if nessary please deal with the NAs.")
    judge_NA_1 <- which(is.na(Data[1,]))
    judge_NA_end <- which(is.na(Data[dim(Data)[1],]))
    if(judge_NA_1){
      for (i in judge_NA_1) {
        jg_n_NA <- min(which(!is.na(Data[,i])))
        Data[1:(jg_n_NA-1),i] <- Data[jg_n_NA,i]
      }
    }

    if(judge_NA_end){
      for (i in judge_NA_end) {
        jg_n_NA <- max(which(!is.na(Data[,i])))
        Data[(jg_n_NA+1):length(Data[,i]),i] <- Data[jg_n_NA,i]
      }
    }
    judge_NA_all <- rle(is.na(Data))
    change_index <- cumsum(judge_NA_all$lengths)
    length_NA <- judge_NA_all$lengths[which(judge_NA_all$values)]
    change_place <- cbind(c(1,change_index)+1, c(change_index, 1))
    change_place <- change_place[-dim(change_place)[1],]
    NA_place <- change_place[which(judge_NA_all$values),]
    Data[which(is.na(Data))] <- do.call(c, apply(cbind(replace_boundary, length_NA), 1, function(re) {seq(re[1], re[2], length.out = 2 + re[3])[-c(1, 2 + re[3])]}))
  }
  dim_Data <- dim(Data)
  judge_rising <- (diff(Data)[-1,] >= -decrease_toleranz)
  return(apply(judge_rising, 2, function(jr) {
    rle_rising <- rle(jr)
    length_rising <- rle_rising$lengths[which(rle_rising$values)]
    length_rising <- length_rising[which(length_rising >= min_rising_duration)]
    return(1 / mean(length_rising))
  }))

}

sig_Peaks_distribution <- function(Data, decrease_toleranz = 0, tile_low = 10, tile_high = 50){
  ## deal with NAs
  if(any(is.na(Data))) interpolat_NA_time(Data)
  dim_Data <- dim(Data)
  diff_Data <- diff(Data)
  peak_location <- rbind(as.matrix(diff(sign(diff_Data))==-2), NA)
  peak_location[!peak_location] <- NA
  Q_peaks <- rbind(NA, as.matrix(Data)) * peak_location
  message("The next message about NAs please irgnore.")
  return(sig_Slope_fdc(Q_peaks, tile_low, tile_high ))
}

mean_Roll_window <- function(Data, roll_size = 7){
  ## deal with NAs
  if(any(is.na(Data))){
    message("The NA in \'Data\' will be linear interpolated with the last and next data, that is not NA, if nessary please deal with the NAs.")
    judge_NA_1 <- which(is.na(Data[1,]))
    judge_NA_end <- which(is.na(Data[dim(Data)[1],]))
    if(judge_NA_1){
      for (i in judge_NA_1) {
        jg_n_NA <- min(which(!is.na(Data[,i])))
        Data[1:(jg_n_NA-1),i] <- Data[jg_n_NA,i]
      }
    }

    if(judge_NA_end){
      for (i in judge_NA_end) {
        jg_n_NA <- max(which(!is.na(Data[,i])))
        Data[(jg_n_NA+1):length(Data[,i]),i] <- Data[jg_n_NA,i]
      }
    }
    judge_NA_all <- rle(is.na(Data))
    change_index <- cumsum(judge_NA_all$lengths)
    length_NA <- judge_NA_all$lengths[which(judge_NA_all$values)]
    change_place <- cbind(c(1,change_index)+1, c(change_index, 1))
    change_place <- change_place[-dim(change_place)[1],]
    NA_place <- change_place[which(judge_NA_all$values),]
    Data[which(is.na(Data))] <- do.call(c, apply(cbind(replace_boundary, length_NA), 1, function(re) {seq(re[1], re[2], length.out = 2 + re[3])[-c(1, 2 + re[3])]}))
  }
  return(rollmean(Data, roll_size))
}

sig_Roll_window <- function(Data, roll_size = 7, max = T){
  Roll_mean <- mean_Roll_window(Data, roll_size)
  ifelse(max, return(apply(Roll_mean, 2, max) * roll_size), return(apply(Roll_mean, 2, min) * roll_size))
}


interpolat_NA_time <- function(Data) {
  message("The NA in \'Data\' will be linear interpolated with the last and next data (in time demision), that is not NA, if nessary please deal with the NAs.")
  judge_NA_1 <- which(is.na(Data[1,]))
  judge_NA_end <- which(is.na(Data[dim(Data)[1],]))
  if(judge_NA_1){
    for (i in judge_NA_1) {
      jg_n_NA <- min(which(!is.na(Data[,i])))
      Data[1:(jg_n_NA-1),i] <- Data[jg_n_NA,i]
    }
  }

  if(judge_NA_end){
    for (i in judge_NA_end) {
      jg_n_NA <- max(which(!is.na(Data[,i])))
      Data[(jg_n_NA+1):length(Data[,i]),i] <- Data[jg_n_NA,i]
    }
  }
  vect_Data <- as.vector(Data)
  judge_NA_all <- rle(is.na(vect_Data))
  change_index <- cumsum(judge_NA_all$lengths)
  length_NA <- judge_NA_all$lengths[which(judge_NA_all$values)]

  change_place <- cbind((append(0,change_index)+1), append(change_index, 1))
  change_place <- change_place[-dim(change_place)[1],]
  NA_place <- matrix(change_place[which(judge_NA_all$values),],,2)
  NA_range <- matrix(0, length(length_NA), 3)
  NA_range[,1] <- vect_Data[NA_place[,1]-1]
  NA_range[,2] <- vect_Data[NA_place[,2]+1]
  NA_range[,3] <- length_NA
  vect_Data[which(is.na(vect_Data))] <- unlist(apply(NA_range, 1, function(rg) seq(rg[1], rg[2], length.out = rg[3] + 2)[-append(1, (rg[3] + 2))]))
  Data[,] <- vect_Data
  return(Data)
}

rev(xts_Qna)
