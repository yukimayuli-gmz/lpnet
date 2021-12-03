#' Only use distance matrix and circular ordering as the input for nnls
#'
#' Using nnls (Non-negative least squares) to calculate weights of splits which are consist with a circular ordering.
#' The input are the distance matrix and the circular ordering.
#' Generally, nnls use \code{A} and \code{b} as input and calculate \code{x} for \code{min(|Ax-b|^2)}.
#' This function will calculate the corresponding matrix \code{A} in Fortran instead of in R which will use less memory and run time.
#' And the vector \code{b} will be get from distance matrix.
#' This function is only a intermediate step of \code{\link{heuristic.method}} and we provide this function in package for test.
#'
#' @param M the distance matrix (the labels of rows and columns are both from 1 to n)
#' @param circular.ordering the circular ordering for a network.
#'
#' @return A list which include:
#'
#' \code{x}, the splits weights estimated by nnls;
#'
#' \code{fitted}, the fitted distance between every taxa pair;
#'
#' \code{deviance}, the residual sum of squares between input distance and fitted distance.
#'
#' (The ordering of splits in result is:
#' \{1|2...n\}, \{1,2|3...n\},..., \{1...n-1|n\}, \{2|1,3...n\},..., \{2...n-1|1,n\},..., \{n-1|1...n-2,n\})
#'
#' @export
nnls.only.use.b <- function(M, circular.ordering) {
  taxa.num<-length(circular.ordering)
  for (i in 1:taxa.num) {
    y<-i
    while (circular.ordering[i]!=i) {
      b0<-M[circular.ordering[i],]
      M[circular.ordering[i],]<-M[y,]
      M[y,]<-b0
      b0<-M[,circular.ordering[i]]
      M[,circular.ordering[i]]<-M[,y]
      M[,y]<-b0
      c0<-circular.ordering[circular.ordering[i]]
      y<-circular.ordering[circular.ordering[i]]<-circular.ordering[i]
      circular.ordering[i]<-c0
    }
  }
  b<-M[row(M)<col(M)]
  MLine <- M <- N <- MDA <- taxa.num*(taxa.num-1)/2
  RNORM <- MODE <- NSETP <- 0
  W <- INDEX <- X <- ZZ <- rep(0, N)
  sol <- .Fortran("nnls_b", Taxanum = as.integer(taxa.num), MLine = as.integer(MLine), B = as.numeric(b),
                  X = as.numeric(X), RNORM = as.numeric(RNORM), W =
                    as.numeric(W), ZZ = as.numeric(ZZ), INDEX =
                    as.integer(INDEX), MODE = as.integer(MODE),
                  NSETP = as.integer(NSETP), PACKAGE="lpnet")
  fitted <- matrix(0,ncol = 1,nrow = M)
  for (j in (2:taxa.num)) {
    for (i in (1:(j-1))) {
      B<-matrix(0,ncol = taxa.num,nrow = taxa.num)
      if(j==taxa.num){
        for (k in (1:i)) {
          B[k,(1:k)]<-0
        }
        for (k in ((i+1):taxa.num)) {
          B[k,(1:i)]<-1
          B[k,((i+1):k)]<-0
        }
      }
      else{
        for (k in (1:i)) {
          B[k,(1:k)]<-0
        }
        for (k in ((i+1):j)) {
          B[k,(1:i)]<-1
          B[k,((i+1):k)]<-0
        }
        for (k in ((j+1):taxa.num)) {
          B[k,(1:i)]<-0
          B[k,((i+1):j)]<-1
          B[k,((j+1):k)]<-0
        }
      }
      fitted[(j-2)*(j-1)/2+i]<-B[col(B)<row(B)] %*% sol$X
    }
  }
  nnls.out <- list(x=sol$X, fitted = fitted, deviance=sol$RNORM^2)
  return(nnls.out)
}
#dyn.load("src/lpnet.dll")
