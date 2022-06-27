#' @title Reshape the time series in water-year
#' @description Reshape the time series in water-year, set the start month of the water year,
#' the number of a water year will be same as the first year-number. And the data before the start month will in the end of the time series.
#' Please make sure the time series begin with the 1. Jan..
#' @param Data [xts]-2D(time, space) Time serious Data in difficult location.
#' @param wateryear_start_m [num]-|2|3|4|5|6|7|8|9|10|11|12| for Feb. ... Dec.
#' @importFrom xts .indexyear .indexyday .indexmon
#' @return By water year reshaped time serious Data [xts]-2D(time, space)
#' @examples
#' ts_water_year(xts_Q, 10)
#' @export
ts_water_year <- function(Data, wateryear_start_m){
  start_y <- .indexyear(Data)[1]
  start_yd <- .indexyday(Data)[1]
  N_preData <- sum(.indexyear(Data) == start_y & .indexmon(Data) < (wateryear_start_m - 1))
  if(start_yd != 0) errorCondition("If you want to use the \"water year\" start with other month, please make sure that the data series begin with 1.Jan..")
  Data[,] <- rbind(as.matrix(Data[-(1:N_preData), ]), as.matrix(Data[(1:N_preData), ]))
  return(Data)
}

#' @title  linear interpolation
#' @description linear interpolated with the last and next data (in time demision)
#' @param Data [xts]-2D(time, space) Time serious Data in difficult location.
#' @return Interpolated time serious Data [xts]-2D(time, space)
#' @examples xts_Q_NA <- xts_Q
#' xts_Q_NA[3:5,1] <- NA
#' rslt <- ts_interpolat_NA_time(xts_Q_NA)
#' @export
ts_interpolat_NA_time <- function(Data) {
  message("The NA in \'Data\' will be linear interpolated with the last and next data (in time demision), that is not NA, if nessary please deal with the NAs.")
  judge_NA_1 <- which(is.na(Data[1,]))
  judge_NA_end <- which(is.na(Data[dim(Data)[1],]))
  if(length(judge_NA_1)){
    for (i in judge_NA_1) {
      jg_n_NA <- min(which(!is.na(Data[,i])))
      Data[1:(jg_n_NA-1),i] <- Data[jg_n_NA,i]
    }
  }

  if(length(judge_NA_end)){
    for (i in judge_NA_end) {
      jg_n_NA <- max(which(!is.na(Data[,i])))
      Data[(jg_n_NA+1):length(Data[,i]),i] <- Data[jg_n_NA,i]
    }
  }
  judge_NA <- which(is.na(Data))
  if(length(judge_NA)){
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
  }
  return(Data)
}

#' @title mean of roll windows
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param roll_size [integer] Roll windows wedith.
#' @importFrom zoo rollmean
#' @return new time series
#' @examples rslt <- ts_Roll_mean(xts_Q)
#' @export
ts_Roll_mean <- function(Data, roll_size = 7){

  return(rollmean(Data, roll_size))
}

#' @title Base flow
#' @description Base flow with digital filter method.
#' @references
#' [Arnold.1995]
#' [eq. 1-2]
#' @param Streamflow [xts]-2D(time, space) Time serious Streamflow in difficult location.
#' @param param_filter [num]-range() Filter parameter that enables the shape of the separation to be altered.
#' @return baseflow [xts]-2D(time, space) Time serious baseflow in difficult location.
#' @examples
#' rslt <- base_Flow.DFM_matrix(xts_Q)
#' @export
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
#' @inheritParams base_Flow.DFM_matrix
#' @param passes_N [integer] The number of pass, see [Arnold.1995].
#' @param negativ_step [logi] TRUE | FALSE for with or without inveser order, see [Arnold.1995].
#' @importFrom xts as.xts
#' @importFrom tidyr replace_na
#' @examples
#' rslt <- base_Flow.DFM(xts_Q)
#' rslt <- base_Flow.DFM(xts_Q, passes_N = 3, negativ_step = FALSE)
#' @export
base_Flow.DFM <- function(Streamflow, param_filter = 0.925, passes_N = 1, negativ_step = TRUE){
  ## deal with NAs
  if(any(is.na(Streamflow))) message("The NA in Streamflow will be replace with 0, if nessary please deal with the NAs.")
  Streamflow[,] <- tidyr::replace_na(as.numeric(Streamflow[,]), 0)
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

fct_filter.DFM <- function(Streamflow, param_filter = 0.925, negativ_step = TRUE){
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

