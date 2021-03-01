#' Read the Taxa block from a Nexus file
#'
#' Read the taxa block from a nexus file, and get the taxa's labels.
#'
#' @param filename a character, the name of the file which will be read for the taxa labels.
#' @param n the number of taxa
#'
#' @return A character vector of the taxa's labels read from the nexus file.
#'
#' @seealso
#' \code{\link{read.nexus.distanceblock}}
#'
#' @export
read.nexus.taxablock<-function(filename,n){
  str<-read.table(filename,sep = ",",fill = TRUE)
  a<-which(str=="BEGIN Taxa;")
  taxaname<-seq(from=0,to=0,length.out = n)
  for (i in 1:n) {
    str2<-str[a+2+i,1]
    str2<-gsub(" ",",",str2,perl = TRUE)
    str2<-strsplit(str2,",")
    str2<-str2[[1]]
    str2<-str2[length(str2)]
    taxaname[i]<-str2
  }
  return(taxaname)
}

#' Read the Distance block from a Nexus file
#'
#' Read the Distance block from a nexus file, and get the symmetry distance matrix.
#' (The distance block should be set as \code{diagonal} and \code{triangle=both}).
#'
#' @param filename a character, the name of the file which will be read for the distance matrix.
#' @param n the number of taxa
#'
#' @return A symmetry distance matrix
#'
#' @seealso
#' \code{\link{read.nexus.taxablock}}
#'
#' @export
read.nexus.distanceblock<-function(filename,n){
  str<-read.table(filename,sep = ",",fill = TRUE)
  a<-which(str=="BEGIN Distances;")
  b<-which(str=="MATRIX")
  a<-b[which(b>a)[1]]
  M<-matrix(0,n,n)
  for (i in 1:n) {
    str2<-str[a+i,1]
    str2<-gsub(" ",",",str2,perl = TRUE)
    str2<-strsplit(str2,",")
    str2<-str2[[1]]
    str2<-str2[(length(str2)-(n-1)):length(str2)]
    str2<-as.numeric(str2)
    M[1:n,i]<-str2
  }
  return(M)
}
