# data_Dresden <- read.csv("Region193_Dresden.csv")
# data_Dresden$date <- as.Date(data_Dresden$date)
# library(xts)
# library(purrr)
# library(plyr)
# library(zoo)
# library(rmatio)
# # library(HMtools)
#
# data29 <- read.mat("33029_daily.mat")
# data29[c(1,6)] <- NULL
# df_data29 <- as.data.frame(data29)
#
# data20 <- read.mat("39020_daily.mat")
# data20[c(1,6)] <- NULL
# df_data20 <- as.data.frame(data20)
#
# df_data29$Date <- seq(as.Date("1999-1-1"), as.Date("2008-12-31"), length.out = 3653)
# xts_Data29 <- xts(df_data29[,1:4], df_data29$Date)
# xts_Data20 <- xts(df_data20, df_data29$Date)
#
# xts_Q <- cbind(xts_Data29$Q, xts_Data20$Q)
# xts_P <- cbind(xts_Data29$P, xts_Data20$P)
# xts_T <- cbind(xts_Data29$T, xts_Data20$T)
# xts_PET <- cbind(xts_Data29$PET, xts_Data20$PET)


# usethis::use_data(xts_Q, xts_P, xts_T, xts_PET)



