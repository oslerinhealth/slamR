# helper functions



#' convert binary codes to decimal integers
#'
#' @param mat A binary matrix the rows of which are the binary codes to be converted.
#' @param LOG Take the logarithm of the decimal integer (default is \code{TRUE});
#' otherwise, set to \code{FALSE}.
#'
#' @return A vector of numbers after conversion
#' @export
#'
#' @examples
#' bin2ind(c(1,1,1))
#' bin2ind(rbind(c(1,0,0),c(1,1,0),c(0,0,1),c(1,0,1),c(0,0,0)))
#' bin2ind(rbind(c(1,0,0),c(1,1,0),c(0,0,1),c(1,0,1),c(0,0,0)),LOG=FALSE)
bin2ind <- function(mat,LOG=FALSE){# This is the log version:
  if(!is.matrix(mat)){mat <- matrix(c(mat),nrow=1)} # convert to a matrix.
  L <- ncol(mat)
  M <- nrow(mat)
  v <- log(as.matrix(2^{(L-1):0},ncol=1))
  diag_v <- matrix(0,nrow=L,ncol=L)
  diag(diag_v) <- v
  tmp_prod <- mat%*%diag_v
  permute_M_vec <- sapply(1:M,function(i){matrixStats::logSumExp(tmp_prod[i,mat[i,]!=0])})
  if (!LOG){
    permute_M_vec <- exp(permute_M_vec)
  }
  #permute_M_vec[rev(order(permute_M_vec))]
  permute_M_vec
}


#' convert decimal integers to binary codes
#'
#' @param x a vector of base 10 integers
#' @return A matrix of binary forms of x (by row)
#' @export
#'
#' @examples
#' decimals <- c(3,5,11,4)
#' m <- binary(decimals)
#' m
binary <- function(x,noBits){
  t(sapply(x,function(number){
    binary_vector = rev(as.numeric(intToBits(number)))
    if(missing(noBits)) {
      return(binary_vector)
    } else {
      binary_vector[-(1:(length(binary_vector) - noBits))]
    }
  }))
}



