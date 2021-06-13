data_Dresden <- read.csv("Region193_Dresden.csv")
data_Dresden$date <- as.Date(data_Dresden$date)
library(xts)
library(purrr)
library(plyr)
library(zoo)
xts_Data <- xts(data_Dresden[1:1000, 1:2], data_Dresden$date[1:1000])

sig_Mean <- function(Data, Time, t_scale = apply.monthly){
  xts_Data <- xts(Data, Time)
  return(t_scale(xts_Data, mean))
}

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

sig_Slope_fdc <- function(Data, tile_low = 33, tile_high = 66){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")

  Quantile_33_66 <- sig_Quantile(Data, c(tile_low, tile_high))
  return((Quantile_33_66[1,] - Quantile_33_66[2,]) / (tile_high/100 - tile_low/100))
}

sig_Runoff_rate <- function(Streamflow, Preciptation){
  return(colSums(Streamflow, na.rm = T) / colSums(Preciptation, na.rm = T))
}

base_Flow.DFM <- function(Streamflow, param_filter = 0.95){
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
  Q_diff <- diff(Streamflow) * (1 + param_filter) * 0.5
  Q_direct <- apply(Q_diff, 2, function(Q_) mat_filter %*% Q_)
  return((Streamflow - Q_direct)[-1, ])
}

sig_Baseflow_index.DFM <- function(Streamflow, param_filter = 0.95){
  Q_base <- base_flow.DFM(Streamflow, param_filter)
  mat_rate <- colSums(Q_base, na.rm = T) / colSums(Streamflow, na.rm = T)
  return(mat_rate)
}

sig_Stream_elas.Sawicz <- function(Streamflow, Preciptation, t_scale = apply.yearly) {
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Preciptation))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  Q_mean_t_scale <- t_scale(Streamflow, mean, na.rm = T)
  P_mean_t_scale <- t_scale(Preciptation, mean, na.rm = T)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  dQ <- diff(Q_mean_t_scale)[-1]
  dP <- diff(P_mean_t_scale)[-1]
  elas_ <- dQ / dP * rep(colMeans(P_mean_t_scale), each = dim_scale[1]) / rep(colMeans(Q_mean_t_scale), dim_scale[1])
  return(elas_)
}

sig_Stream_elas.Sanka <- function(Streamflow, Preciptation, t_scale = apply.yearly) {
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Preciptation))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  Q_mean_t_scale <- t_scale(Streamflow, mean, na.rm = T)
  P_mean_t_scale <- t_scale(Preciptation, mean, na.rm = T)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  Q_mean <- rep(colMeans(Q_mean_t_scale), dim_scale[1])
  P_mean <- rep(colMeans(P_mean_t_scale), each = dim_scale[1])
  elas_ <- (Q_mean_t_scale - Q_mean) / (P_mean_t_scale - P_mean) * P_mean / Q_mean
  return(elas_[-1,])
}

sig_Frequency_threshold <- function(Data, higher = T, threshold = c(1, 3, 9), times_median = T){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")

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

sig_Duration_threshold <- function(Data, higher = T, threshold = c(1, 3, 9), times_median = T, dur_type = mean){
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

sig_Autocor <- function(Data){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  time_N <- dim(Data)[1]
  return(apply(Data, 2, function(D) cor(D[-1], D[-time_N], na.rm = T)))
}

sig_Peaks_distribution <- function(Data, decrease_toleranz = 0, tile_low = 10, tile_high = 50){
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
