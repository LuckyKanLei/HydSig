#' @title Mean
#' @description  Mean of data, e.g. streamflow, precipitation, temperature and so on.
#' @param Data [xts]-2D(time, space) Time serious Data in difficult location.
#' @param t_scale [function / NULL] |apply.daily| for daily, |apply.weekly| for weekly, |apply.monthly| for monthly,
#' |apply.quarterly| for quarterly, |apply.yearly| for yearly, |NULL| for all of the time.
#' @importFrom xts apply.daily apply.weekly apply.monthly apply.quarterly apply.yearly
#' @return [xts / array]-2D(time, space) Mean of data.
#' @examples
#' rslt <- sig_Mean(xts_Q)
#' rslt <- sig_Mean(xts_Q, t_scale = NULL)
#' @export
sig_Mean <- function(Data, t_scale = apply.monthly){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  if(is.null(t_scale)) return(colMeans(Data, na.rm = TRUE))
  else return(t_scale(Data, colMeans, na.rm = TRUE))
}

#' @title  Quantile / Percentile decreasing
#' @description Quantile / Percentile of data in decreasing, e.g. streamflow, precipitation, temperature and so on.
#' The Q5 is the high flow, Q95 is the low flow.
#' For the other scale (if you want) can use apply.yearly(Data, sig_Quantile, tile = c(5, 95)) or apply.weekly apply.monthly apply.quarterly apply.yearly.
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param tile [num / vector]-range(0-100) The tiles (positions).
#' @return [xts / array]-2D(tile, space) The data in the tiles.
#' @examples
#' rslt <- sig_Quantile(xts_Q)
#' rslt <- sig_Quantile(xts_Q, tile = c(10, 50, 90))
#' @export
sig_Quantile <- function(Data, tile = 95){
  if(any(tile < 0 | tile > 100)) errorCondition("Please give the tile between 0 and 100.")
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  tile <- tile / 100.
  tile_N <- length(tile)
  sort_Data <- apply(Data, 2, sort, decreasing = TRUE, na.last = TRUE)
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
#' @examples rslt <- sig_Runoff_ratio(xts_Q, xts_P)
#' @export
sig_Runoff_ratio <- function(Streamflow, Preciptation){
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(c(Streamflow, Preciptation)))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  return(colSums(Streamflow, na.rm = TRUE) / colSums(Preciptation, na.rm = TRUE))
}


#' @title Baseflow index
#' @inheritParams base_Flow.DFM
#' @param t_scale [function / NULL] |apply.daily| for daily, |apply.weekly| for weekly, |apply.monthly| for monthly,
#' |apply.quarterly| for quarterly, |apply.yearly| for yearly, |NULL| for all of the time.
#' @return Baseflow index [xts]-2D(1, space) Time serious baseflow in difficult location.
#' @examples sig_Baseflow_index.DFM(xts_Q, negativ_step = FALSE, passes_N = 3, t_scale = NULL)
#' sig_Baseflow_index.DFM(xts_Q)
#' @export
sig_Baseflow_index.DFM <- function(Streamflow, param_filter = 0.925, passes_N = 1, negativ_step = TRUE, t_scale = apply.yearly){
  if(any(is.na(Streamflow))) message("The NA in Streamflow will be replace with 0, if nessary please deal with the NAs.")
  Streamflow[,] <- tidyr::replace_na(as.numeric(Streamflow[,]), 0)
  Q_base <- base_Flow.DFM(Streamflow, param_filter, passes_N, negativ_step)
  if(is.null(t_scale)) mat_rate <- colSums(Q_base, na.rm = TRUE) / colSums(Streamflow, na.rm = TRUE)
  else {
    BQ <- cbind(B = Q_base, Q = Streamflow)
    n_gauge <- dim(BQ)[2] / 2
    mat_temp <- t_scale(BQ, colMeans, na.rm = TRUE)
    mat_rate <- mat_temp[,1:n_gauge] / mat_temp[, (n_gauge + 1):(n_gauge * 2)]
  }
  return(mat_rate)
}


#' @title Streamflow-precipitation elasticity
#' @param Streamflow [xts]-2D(time, space) Time serious Streamflow in difficult location.
#' @param Preciptation [xts]-2D(time, space) Time serious Preciptation in difficult location.
#' @param t_scale [function] |apply.quarterly| for quarterly, |apply.yearly| for yearly.
#' @param wateryear_start_m [num]-|2|3|4|5|6|7|8|9|10|11|12| for Feb. ... Dec.
#' @param fct_stat [function] statistic function like mean, media, max and so on
#' @return stream elasticity [xts]-2D(time, space) in time scale
#' @examples rslt <- sig_Stream_elas.Sawicz(xts_Q, xts_P)
#' rslt <- median(sig_Stream_elas.Sawicz(xts_Q, xts_P, wateryear_start_m = 10))
#' @export
sig_Stream_elas.Sawicz <- function(Streamflow, Preciptation, t_scale = apply.yearly, wateryear_start_m = 1, fct_stat = mean) {
  if(wateryear_start_m != 1 && as.character(substitute(t_scale)) != "apply.yearly") message("Water year is only for yearly availability.")
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Preciptation))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  if(wateryear_start_m > 1) {
    Streamflow <- ts_water_year(Streamflow, wateryear_start_m)
    Preciptation <- ts_water_year(Preciptation, wateryear_start_m)
  }
  Q_mean_t_scale <- t_scale(Streamflow, colMeans, na.rm = TRUE)
  P_mean_t_scale <- t_scale(Preciptation, colMeans, na.rm = TRUE)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  dQ <- diff(Q_mean_t_scale)
  dP <- diff(P_mean_t_scale)
  elas_ <- (dQ / dP * rep(colMeans(P_mean_t_scale), each = dim_scale[1]) / rep(colMeans(Q_mean_t_scale), dim_scale[1]))[-1,]
  return(apply(elas_, 2, fct_stat, na.rm = T))
}


#' @rdname sig_Stream_elas.Sawicz
#' @examples rslt <- sig_Stream_elas.Sanka(xts_Q, xts_P)
#' rslt <- median(sig_Stream_elas.Sanka(xts_Q, xts_P, wateryear_start_m = 10))
#' @export
sig_Stream_elas.Sanka <- function(Streamflow, Preciptation, t_scale = apply.yearly, wateryear_start_m = 1, fct_stat = mean) {
  if(wateryear_start_m != 1 && as.character(substitute(t_scale)) != "apply.yearly") message("Water year is only for yearly availability.")
  dim_Q <- dim(Streamflow)
  dim_P <- dim(Preciptation)
  if(dim_Q[1] != dim_P[1] | dim_Q[2] != dim_P[2]) errorCondition("Please make sure Streamflow and Preciptation have the same dimsion.")
  if(any(is.na(Streamflow)) | any(is.na(Preciptation))) message("The NA in Streamflow and Preciptation will be droped, if nessary please deal with the NAs.")
  Q_mean_t_scale <- t_scale(Streamflow, colMeans, na.rm = TRUE)
  P_mean_t_scale <- t_scale(Preciptation, colMeans, na.rm = TRUE)
  dim_scale <- dim(Q_mean_t_scale)
  if(dim_scale[1] < 2) errorCondition("The average Streamflow and Preciptation in the number of time intervall must bigger than 2.")
  Q_mean <- rep(colMeans(Streamflow, na.rm = TRUE), each = dim_scale[1])
  P_mean <- rep(colMeans(Preciptation, na.rm = TRUE), each = dim_scale[1])
  elas_ <- (Q_mean_t_scale - Q_mean) / (P_mean_t_scale - P_mean) * P_mean / Q_mean
  return(apply(elas_, 2, fct_stat, na.rm = T))
}

#' @title Frequency or Duration by specified threshold
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param higher [logi] TRUE | FALSE for the value bigger or smaller than the threshold.
#' @param threshold [num / vector] one value or a vactor of threshold, If times_mean = TRUE, it should be the times with mean value;
#' if not it can be the concrete value in the range of \'Data\' and without times the mean value.
#' @param times_mean [logi] TRUE | FALSE for use or not use the times with mean value.
#' @return frequency [num]
#' @examples  rslt <- sig_Frequency_threshold(xts_Q)
#' rslt <- sig_Frequency_threshold(xts_Q, higher = FALSE, threshold = c(0.1, 1), times_mean = FALSE)
#' @export
sig_Frequency_threshold <- function(Data, higher = TRUE, threshold = 1, times_mean = TRUE){
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
      if(dim(threshold)[2] != space_N) errorCondition("The second dimesion of \'threshold\' must same as the \'Data\'.")
      else threshold_N <- dim(threshold)[1]}
  }

  dim_Data <- dim(Data)
  mat_freq <- matrix(0, threshold_N, space_N)
  colnames(mat_freq) <- colnames(Data)
  if(higher){
    for (i in 1:threshold_N) {
      mat_Diff <- Data - rep(threshold[i,], each = time_N)
      mat_freq[i,] <- colSums(mat_Diff > 0, na.rm = TRUE) / colSums(!is.na(Data), na.rm = TRUE)
    }
  } else {
    for (i in 1:threshold_N) {
      mat_Diff <- Data - rep(threshold[i,], each = time_N)
      mat_freq[i,] <- colSums(mat_Diff <= 0, na.rm = TRUE) / colSums(!is.na(Data), na.rm = TRUE)
    }
  }
  return(mat_freq)
}

#' @describeIn sig_Frequency_threshold Frequency or Duration by specified threshold
#' @importFrom tidyr replace_na
#' @param dur_type [function] mean | max | min for the stattistical methods for the qualifiede Durations.
#' @examples  rslt <- sig_Duration_threshold(xts_Q)
#' rslt <- sig_Duration_threshold(xts_Q, higher = FALSE, threshold = c(0.1, 1),
#' times_mean = FALSE, dur_type = max)
#' @export
sig_Duration_threshold <- function(Data, higher = TRUE, threshold = 1, times_mean = TRUE, dur_type = mean){
  ## deal with NAs
  Data[,] <- tidyr::replace_na(as.numeric(Data[,]), ifelse(higher, min(Data), max(Data)))
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

#' @title  Snowday ratio
#' @description the ration of the day with snow, with the Temperatur methods
#' @param Preciptation [xts]-2D(time, space) Time serious Preciptation in difficult location.
#' @param Temperatur [xts]-2D(time, space) Time serious Temperatur in difficult location.
#' @param snow_threshold Temperatur threshold for the snow day
#' @examples rslt <- sig_Snowday_ratio(xts_P, xts_Q)
#' @export
sig_Snowday_ratio <- function(Preciptation, Temperatur, snow_threshold = 2){
  dim_T <- dim(Temperatur)
  dim_P <- dim(Preciptation)
  if(dim_T[1] != dim_P[1] | dim_T[2] != dim_P[2]) errorCondition("Please make sure Temperatur and Preciptation have the same dimsion.")
  if(any(is.na(Preciptation)) | any(is.na(Temperatur))) message("The NA in Temperatur and Preciptation will be droped, if nessary please deal with the NAs.")

  judge_P <- Preciptation > 0
  judge_S <- judge_P & (Temperatur < snow_threshold)
  S_ration <- colSums(judge_S, na.rm = TRUE) / colSums(judge_P, na.rm = TRUE)
  S_ration[is.na(S_ration)] <- 0
  return(S_ration)
}

#' @title Autocorrelation
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param time_lag [integer] The number of lag.
#' @return Autocorrelation [array]-2D(1, space)
#' @examples rslt <- sig_Autocor(xts_Q)
#' @export
sig_Autocor <- function(Data, time_lag = 1){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  time_N <- dim(Data)[1]
  return(apply(Data, 2, function(D) {
    D <- as.numeric(D[which(!is.na(D))])
    cor(D[-((length(D) - time_lag + 1):length(D))], D[-(1:time_lag)])}))
}


#' @title Rising limb density
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param decrease_toleranz [num] Toleranz for the small decrease.
#' @param min_rising_duration [integer] The smallest duration for a decrease.
#' @importFrom stats cor
#' @importFrom tidyr replace_na
#' @return sig_Rising_limb_density [num] [1/step]
#' @examples rslt <- sig_Rising_limb_density(xts_Q)
#' @export
sig_Rising_limb_density <- function(Data, decrease_toleranz = 0, min_rising_duration = 1){
  ## deal with NAs
  if(any(is.na(Data))) message("The NA in Data will be replace with 0, if nessary please deal with the NAs.")
  Data[,] <- tidyr::replace_na(as.numeric(Data[,]), 0)
  dim_Data <- dim(Data)
  judge_rising <- (diff(Data)[-1,] >= -decrease_toleranz)
  return(apply(judge_rising, 2, function(jr) {
    rle_rising <- rle(jr)
    length_rising <- rle_rising$lengths[which(rle_rising$values)]
    length_rising <- length_rising[which(length_rising >= min_rising_duration)]
    return(1 / mean(length_rising))
  }))

}



#' @title Slope of the flow duration curve
#' @description Slope of the flow duration curve of the data mostly streamflow.
#' For the other scale can use apply.yearly(Data, sig_Quantile, tile = c(5, 95)) or apply.weekly apply.monthly apply.quarterly apply.yearly.
#' @references [Sawicz.2011]
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param tile_low [num]-range(0, 100) The lower tile, but the data is the bigger, because the FDC is in decreasing.
#' @param tile_high [num]-range(0, 100) The higher tile, but the data is the smaller, because the FDC is in decreasing.
#' @param log_fit |T| for mit ln(Q), |F| for ohne ln(Q)
#' @return [array]-2D(1, space) The data in the tiles. Or [xts]-2D(sum_time, space) for other scale with apply.XXXly.
#' @examples rslt <- sig_Slope_fdc(xts_Q)
#' rslt <- sig_Slope_fdc(xts_Q, 5, 50, FALSE)
#' rslt <- xts::apply.yearly(xts_Q, sig_Slope_fdc)
#' @export
sig_Slope_fdc <- function(Data, tile_low = 33, tile_high = 66, log_fit = T){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")

  Quantile_33_66 <- sig_Quantile(Data, c(tile_low, tile_high))
  if(log_fit) return((log(Quantile_33_66[1,]) - log(Quantile_33_66[2,])) / (tile_low/100 - tile_high/100))
  else return((Quantile_33_66[1,] - Quantile_33_66[2,]) / (tile_low/100 - tile_high/100))
}

#' @title Slope of distribution of peaks
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param decrease_toleranz [num] Toleranz for the small decrease.
#' @param tile_low,tile_high [num]-range(0-100) The title range for slope of the curve.
#' @return sig_Peaks_distribution [num] [1/step]
#' @examples rslt <- sig_Peaks_distribution(xts_Q)
#' @export
sig_Peaks_distribution <- function(Data, decrease_toleranz = 0, tile_low = 10, tile_high = 50){
  ## deal with NAs
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")
  Data <- replace(Data, is.na(Data), 0)

  dim_Data <- dim(Data)
  diff_Data <- diff(Data)
  peak_location <- rbind(as.matrix(diff(sign(diff_Data))==-2), NA)
  peak_location[!peak_location] <- NA
  Q_peaks <- rbind(NA, as.matrix(Data)) * peak_location
  message("The next message about NAs please irgnore.")
  return(sig_Slope_fdc(Q_peaks, tile_low, tile_high ))
}

