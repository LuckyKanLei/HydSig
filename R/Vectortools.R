#' @title Finding the minimum value of the corresponding position of two equal large vectors
#' @param VctA vector, numic, first Vector
#' @param VctB vector, numic, second Vector
#' @return vector, result
#' @examples
#' A = c(1, 2, 3)
#' B = c(2, 2, 2)
#' minVector(A, B)
#' @export
minVector <- function(VctA, VctB){
  VctC = VctA - VctB
  VctD = VctB
  VctD[which(VctC < 0.0)] = VctA[which(VctC < 0.0)]
  return(VctD)
}


#' @rdname minVector
#' @examples
#' maxVector(A, B)
#' @export
maxVector <- function(VctA, VctB){
  VctC = VctA - VctB
  VctD = VctA
  VctD[which(VctC < 0.0)] = VctB[which(VctC < 0.0)]
  return(VctD)
}

#' @title  Finding the minimum / maximum value of one vectors and one wert
#' @param WrtA numic
#' @param VctB vector, numic
#' @return vector, result
#' @examples
#' a = 2
#' B = c(1, 2, 3)
#' maxVector(a, B)
#' @export
minSVector <- function(WrtA, VctB){
  VctC = WrtA - VctB
  VctD = VctB
  VctD[which(VctC < 0.0)] = WrtA
  return(VctD)
}

#' @rdname minSVector
#' @examples
#' maxVector(a, B)
#' @export
maxSVector <- function(WrtA, VctB){
  VctC = WrtA - VctB
  VctD = VctB
  VctD[which(VctC > 0.0)] = WrtA
  return(VctD)
}


