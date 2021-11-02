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
  if(any(is.na(Data))) interpolat_NA_time(Data)
  dim_Data <- dim(Data)
  diff_Data <- diff(Data)
  peak_location <- rbind(as.matrix(diff(sign(diff_Data))==-2), NA)
  peak_location[!peak_location] <- NA
  Q_peaks <- rbind(NA, as.matrix(Data)) * peak_location
  message("The next message about NAs please irgnore.")
  return(sig_Slope_fdc(Q_peaks, tile_low, tile_high ))
}
