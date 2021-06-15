#' @title Slope of the flow duration curve
#' @description Slope of the flow duration curve of the data mostly streamflow.
#' For the other scale can use apply.yearly(Data, sig_Quantile, tile = c(5, 95)) or apply.weekly apply.monthly apply.quarterly apply.yearly.
#' @references [Sawicz.2011]
#' @param Data [xts / array]-2D(time, space) Time serious Data in difficult location.
#' @param tile_low [num]-range(0, 100) The lower tile, but the data is the bigger, because the FDC is in decreasing.
#' @param tile_high [num]-range(0, 100) The higher tile, but the data is the smaller, because the FDC is in decreasing.
#' @param log_fit |T| for mit ln(Q), |F| for ohne ln(Q)
#' @return [array]-2D(1, space) The data in the tiles. Or [xts]-2D(sum_time, space) for other scale with apply.XXXly.
#' @examples sig_Slope_fdc(xts_Q)
#' sig_Slope_fdc(xts_Q, 5, 50, F)
#' apply.yearly(xts_Q, sig_Slope_fdc)
sig_Slope_fdc <- function(Data, tile_low = 33, tile_high = 66, log_fit = T){
  if(any(is.na(Data))) message("The NA in \'Data\' will be droped, if nessary please deal with the NAs.")

  Quantile_33_66 <- sig_Quantile(Data, c(tile_low, tile_high))
  if(log_fit) return((log(Quantile_33_66[1,]) - log(Quantile_33_66[2,])) / (tile_low/100 - tile_high/100))
  else return((Quantile_33_66[1,] - Quantile_33_66[2,]) / (tile_low/100 - tile_high/100))
}
