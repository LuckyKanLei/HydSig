data_Dresden <- read.csv("Region193_Dresden.csv")
data_Dresden$date <- as.Date(data_Dresden$date)
library(xts)

sig_mean <- function(Data, Time, f_type = apply.monthly){
  xts_Data <- xts(Data, Time)
  return(f_type(xts_Data, mean))
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
return((smaller_Quan + bigger_Quan) / 2)
}
sig_slope_fdc(data_[,1:2])
sig_slope_fdc <- function(Data, tile1 = 33, tile2 = 66){
  Quantile_33_66 <- sig_quantile(Data, c(tile1, tile2))
  return((Quantile_33_66[2,] - Quantile_33_66[1,]) / (tile2 - tile1))
}
