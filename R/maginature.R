data_Dresden <- read.csv("Region193_Dresden.csv")
data_Dresden$date <- as.Date(data_Dresden$date)
library(xts)
library(purrr)
xts_Data <- xts(data_Dresden[1:1000, 1:2], data_Dresden$date[1:1000])

sig_mean <- function(Data, Time, t_scale = apply.monthly){
  xts_Data <- xts(Data, Time)
  return(t_scale(xts_Data, mean))
}

sig_quantile <- function(Data, tile = c(5, 95)){
  if(any(tile < 0 | tile > 100)) errorCondition("Please give the tile between 0 and 100.")
  tile <- tile / 100.
  tile_N <- length(tile)
  sort_Data <- apply(Data, 2, sort, na.last = T)
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

sig_slope_fdc <- function(Data, tile_low = 33, tile_high = 66){
  Quantile_33_66 <- sig_quantile(Data, c(tile_low, tile_high))
  return((Quantile_33_66[2,] - Quantile_33_66[1,]) / (tile_high - tile_low))
}

sig_runoff_rate <- function(Streamflow, Preciptation){
  return(colSums(Streamflow, na.rm = T) / colSums(Preciptation, na.rm = T))
}

base_flow.DFM <- function(Streamflow, param_filter = 0.95){
  dim_Q <- dim(Streamflow)
  time_N <- dim_Q[1]
  mat_filter <- matrix(param_filter^(rep(1:time_N,time_N) - rep(1:time_N, each=time_N)), time_N)
  mat_filter[upper.tri(mat_filter)] <- 0
  if(any(is.na(Streamflow))) {
    message("The Null (NA) values in Streamflow will replace by 0, if nessary please deal with the NAs. And the results in the same place will in NA seted.")
    Streamflow0 <- replace(Streamflow, is.na(Streamflow), 0)
    Q_diff <- (Streamflow0 - as.matrix(rbind(Streamflow[1, ], Streamflow0[1:(dim_Q[1]-1), ]))) * (1 + param_filter) * 0.5
    Q_direct <- apply(Q_diff, 2, function(Q_) mat_filter %*% Q_)
    return((Streamflow - Q_direct)[-1, ])
  }
  Q_diff <- (Streamflow - as.matrix(rbind(Streamflow[1, ], Streamflow[1:(dim_Q[1]-1), ]))) * (1 + param_filter) * 0.5
  Q_direct <- apply(Q_diff, 2, function(Q_) mat_filter %*% Q_)
  return((Streamflow - Q_direct)[-1, ])
}

sig_baseflow_index.DFM <- function(Streamflow, param_filter = 0.95){
  Q_base <- base_flow.DFM(Streamflow, param_filter)
  mat_rate <- colSums(Q_base, na.rm = T) / colSums(Streamflow, na.rm = T)
  return(mat_rate)
}

sig_stream_elas.Sawicz <- function(Streamflow, Preciptation, t_scale = apply.yearly) {
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Streamflow))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  Q_mean_t_scale <- t_scale(Streamflow, mean, na.rm = T)
  P_mean_t_scale <- t_scale(Preciptation, mean, na.rm = T)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  dQ <- Q_mean_t_scale - as.matrix(rbind(Q_mean_t_scale[1, ], Q_mean_t_scale[1:(dim_scale[1]-1), ]))
  dP <- P_mean_t_scale - as.matrix(rbind(P_mean_t_scale[1, ], P_mean_t_scale[1:(dim_scale[1]-1), ]))
  elas_ <- dQ / dP * rep(colMeans(P_mean_t_scale), each = dim_scale[1]) / rep(colMeans(Q_mean_t_scale), dim_scale[1])
  return(elas_[-1,])
}

sig_stream_elas.Sanka <- function(Streamflow, Preciptation, t_scale = apply.yearly) {
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Streamflow))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  Q_mean_t_scale <- t_scale(Streamflow, mean, na.rm = T)
  P_mean_t_scale <- t_scale(Preciptation, mean, na.rm = T)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  Q_mean <- rep(colMeans(Q_mean_t_scale), dim_scale[1])
  P_mean <- rep(colMeans(P_mean_t_scale), each = dim_scale[1])
  elas_ <- (Q_mean_t_scale - Q_mean) / (P_mean_t_scale - P_mean) * P_mean / Q_mean
  return(elas_[-1,])
}

sig_frequency_threshold <- function(Data, higher = T, threshold = c(1, 3, 9), times_median = T){
  time_N <- dim(Data)[1]
  space_N <- dim(Data)[2]
  if(times_median) {
    threshold_N <- length(threshold)
    threshold <- matrix(rep(apply(Data, 2, median), each = threshold_N) * rep(threshold, space_N), threshold_N, space_N)
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


sig_duration_threshold <- function(Data, higher = T, threshold = c(1, 3, 9), times_median = T, dur_type = mean){
  time_N <- dim(Data)[1]
  space_N <- dim(Data)[2]
  if(times_median) {
    threshold_N <- length(threshold)
    threshold <- matrix(rep(apply(Data, 2, median), each = threshold_N) * rep(threshold, space_N), threshold_N, space_N)
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


xts_Data
a <- sig_stream_elas.Sanka(xts_Data, xts_Data+10)
sig_duration_threshold(xts_Data)
apply(a, 2, median)
str(xts_Data$q_obs)
colm

str(rle(as.vector(xts_Data$q_obs[1:100]) > 139))
