#' @title Reshape the time series in water-year
#' @description Reshape the time series in water-year, set the start month of the water year,
#' the number of a water year will be same as the first year-number. And the data before the start month will in the end of the time series.
#' Please make sure the time series begin with the 1. Jan..
#' @param Data [xts]-2D(time, space) Time serious Data in difficult location.
#' @param wateryear_start_m [num]-|2|3|4|5|6|7|8|9|10|11|12| for Feb. ... Dec.
#' @importFrom xts .indexyear .indexyday .indexmon
#' @return By water year reshaped time serious Data [xts]-2D(time, space)
#' @examples
#' water_year(xts_Q, 10)
#' @export
water_year <- function(Data, wateryear_start_m){
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
#' rslt <- interpolat_NA_time(xts_Q_NA)
#' @export
interpolat_NA_time <- function(Data) {
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


