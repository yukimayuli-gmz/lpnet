#calculate sum of quartets weights for 4 taxa set
sum_index<-function(a,b,c,d,M){
  la<-length(a)
  lb<-length(b)
  lc<-length(c)
  ld<-length(d)
  A<-array(0,dim = c(la,lb,lc,ld))
  for (i in (1:la)) {
    for (j in (1:lb)) {
      for (k in (1:lc)) {
        for (l in (1:ld)) {
          #ac+bd-ad-bc
          A[i,j,k,l]<-0.5*(M[a[i],c[k]]+M[b[j],d[l]]-M[a[i],d[l]]-M[b[j],c[k]])
        }
      }
    }
  }
  s<-sum(A)
  return(s)
}
#calculate the complement for set a in n
complement_c<-function(n,a){
  z<-c()
  for (i in (1:n)) {
    if(!(i %in% a)){
      z<-c(z,i)
    }
  }
  return(z)
}
#calculate objective coefficients for lp
make_lp_c<-function(n,M,tredge,taxa){
  s<-0
  lp_c<-c()
  for (i in ((n+1):(2*n-3))) {
    for (j in ((i+1):(2*n-2))) {
      s<-s+1
      if(i==(n+1)){
        t1<-tredge[tredge[,1]==i][4]
        t2<-tredge[tredge[,1]==i][5]
        t3<-tredge[tredge[,1]==i][6]
        if(taxa[[j]][1] %in% taxa[[t1]]){
          a<-taxa[[t2]]
          b<-taxa[[t3]]
        }
        if(taxa[[j]][1] %in% taxa[[t2]]){
          a<-taxa[[t3]]
          b<-taxa[[t1]]
        }
        if(taxa[[j]][1] %in% taxa[[t3]]){
          a<-taxa[[t1]]
          b<-taxa[[t2]]
        }
        c<-taxa[[tredge[tredge[,1]==j][3]]]
        d<-taxa[[tredge[tredge[,1]==j][4]]]
      }
      if(i!=(n+1)){
        s1<-tredge[tredge[,1]==i][3]
        s2<-tredge[tredge[,1]==i][4]
        if(taxa[[j]][1] %in% taxa[[s1]]){
          a<-taxa[[s2]]
          b<-complement_c(n,c(taxa[[s1]],taxa[[s2]]))
          c<-taxa[[tredge[tredge[,1]==j][3]]]
          d<-taxa[[tredge[tredge[,1]==j][4]]]
        }
        if(taxa[[j]][1] %in% taxa[[s2]]){
          a<-complement_c(n,c(taxa[[s1]],taxa[[s2]]))
          b<-taxa[[s1]]
          c<-taxa[[tredge[tredge[,1]==j][3]]]
          d<-taxa[[tredge[tredge[,1]==j][4]]]
        }
        if((!(taxa[[j]][1] %in% taxa[[s1]]))&(!(taxa[[j]][1] %in% taxa[[s2]]))){
          a<-taxa[[tredge[tredge[,1]==i][3]]]
          b<-taxa[[tredge[tredge[,1]==i][4]]]
          c<-taxa[[tredge[tredge[,1]==j][3]]]
          d<-taxa[[tredge[tredge[,1]==j][4]]]
        }
      }
      lp_c[s]<-sum_index(a,b,d,c,M)-sum_index(a,b,c,d,M)
    }
  }
  return(lp_c)
}
#calculate the constraint coefficients for lp
make_lp_matrix<-function(n){
  s<-0
  lp_matrix<-matrix(0,nrow = 4*(n-2)*(n-3)*(n-4)/6 , ncol = (n-2)*(n-3)/2)
  for (i in ((n+1):(2*n-4))) {
    for (j in ((i+1):(2*n-3))) {
      for (k in ((j+1):(2*n-2))) {
        lp_matrix[(4*s+1),((i-n-1)*(n-3)+j-n-1-(i-n-1)*(i-n)/2)]<- 1
        lp_matrix[(4*s+1),((i-n-1)*(n-3)+k-n-1-(i-n-1)*(i-n)/2)]<- 1
        lp_matrix[(4*s+1),((j-n-1)*(n-3)+k-n-1-(j-n-1)*(j-n)/2)]<- 1
        lp_matrix[(4*s+2),((i-n-1)*(n-3)+j-n-1-(i-n-1)*(i-n)/2)]<- -1
        lp_matrix[(4*s+2),((i-n-1)*(n-3)+k-n-1-(i-n-1)*(i-n)/2)]<- -1
        lp_matrix[(4*s+2),((j-n-1)*(n-3)+k-n-1-(j-n-1)*(j-n)/2)]<- 1
        lp_matrix[(4*s+3),((i-n-1)*(n-3)+j-n-1-(i-n-1)*(i-n)/2)]<- -1
        lp_matrix[(4*s+3),((i-n-1)*(n-3)+k-n-1-(i-n-1)*(i-n)/2)]<- 1
        lp_matrix[(4*s+3),((j-n-1)*(n-3)+k-n-1-(j-n-1)*(j-n)/2)]<- -1
        lp_matrix[(4*s+4),((i-n-1)*(n-3)+j-n-1-(i-n-1)*(i-n)/2)]<- 1
        lp_matrix[(4*s+4),((i-n-1)*(n-3)+k-n-1-(i-n-1)*(i-n)/2)]<- -1
        lp_matrix[(4*s+4),((j-n-1)*(n-3)+k-n-1-(j-n-1)*(j-n)/2)]<- -1
        s<-s+1
      }
    }
  }
  return(lp_matrix)
}
#calculate the right hand side of the constraints
make_lp_rhs<-function(n){
  lp_rhs<-c()
  for (i in (1:((n-2)*(n-3)*(n-4)/6))) {
    lp_rhs[4*i-3]<-2
    lp_rhs[4*i-2]<-0
    lp_rhs[4*i-1]<-0
    lp_rhs[4*i]<-0
  }
  return(lp_rhs)
}
#give the directions of the constraints
make_lp_dir<-function(n){
  lp_dir<-c()
  for (i in (1:(4*(n-2)*(n-3)*(n-4)/6))) {
    lp_dir[i]<-"<="
  }
  return(lp_dir)
}
#calculate the position of taxa set b in ordering a
which_ordering_in_taxa<-function(a,b){
  x<-c()
  t<-c()
  l<-length(b)
  for (i in 1:l) {
    x[i]<-which(a==b[i])
  }
  t[1]<-min(x)
  t[2]<-max(x)
  return(t)
}
#flip the part of ordering a which position is t
reverse<-function(a,t){
  a[t[2]:t[1]]<-a[t[1]:t[2]]
  return(a)
}
#calculate the ordering of lpnet
lp_opt<-function(ordering,edge,tredge,taxa){
  n<-length(ordering)
  for (i in ((n+2):(2*n-2))) {
    t<-which_ordering_in_taxa(ordering,taxa[[i]])
    i0<-tredge[tredge[,2]==i][1]
    f<-edge[((i0-n-1)*(n-3)+i-n-1-(i0-n-1)*(i0-n)/2)]
    if(f>0.1){
      ordering<-reverse(ordering,t)
    }
  }
  return(ordering)
}

#make the corresponding matrix for nnls
make.matrix <- function(n){
  A <- matrix(nrow = (n*(n-1)/2),ncol = (n*(n-1)/2))
  a<-0
  for (i in (1:(n-1))) {
    for (j in ((i+1):n)) {
      a<-a+1
      B <- matrix(nrow = n,ncol = n)
      if(i==1){
        for (k in (1:(j-1))) {
          B[k,(1:k)]<-0
        }
        for (k in (j:n)) {
          B[k,(1:(j-1))]<-1
          B[k,(j:k)]<-0
        }
      }
      else{
        for (k in (1:(i-1))) {
          B[k,(1:k)]<-0
        }
        for (k in (i:(j-1))) {
          B[k,(1:(i-1))]<-1
          B[k,(i:k)]<-0
        }
        for (k in (j:n)) {
          B[k,(1:(i-1))]<-0
          B[k,(i:(j-1))]<-1
          B[k,(j:k)]<-0
        }
      }
      s<-1
      for (k in (2:n)) {
        for (l in (1:(k-1))) {
          A[s,a] <- B[k,l]
          s<-s+1
        }
      }
    }
  }
  return(A)
}
#write the content of nexus file
make.nexus_split_block_order<-function(o,p,k,taxaname=NULL){
  a_residue<-nnls(make.matrix(k),p)
  a<-a_residue$x
  la<-length(which(a>0))
  x<-"#nexus"
  x<-c(x,"BEGIN Taxa;")
  x<-c(x,paste("DIMENSIONS",gsub(" ","",paste("ntax=",k,";"),perl = TRUE)))
  x<-c(x,"TAXLABELS")
  for (i in (1:k)) {
    if(is.null(taxaname[o[i]])){
      x<-c(x,paste(gsub(" ","",paste("[",i,"]"),perl = TRUE),sub(" ","",paste("taxa",o[i]),perl = TRUE)))
    }else{
      x<-c(x,paste(gsub(" ","",paste("[",i,"]"),perl = TRUE),gsub(" ","",paste("'",taxaname[o[i]],"'"),perl = TRUE)))
    }
  }
  x<-c(x,";")
  x<-c(x,"END; [Taxa]")
  x<-c(x,"BEGIN Splits;")
  x<-c(x,paste("DIMENSIONS",gsub(" ","",paste("ntax=",k),perl = TRUE),gsub(" ","",paste("nsplits=",la,";"),perl = TRUE)))
  x<-c(x,"FORMAT labels=no weights=yes confidences=no intervals=no;")
  x<-c(x,"MATRIX")
  s<-0
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      s<-s+1
      if(a[s]>0){
        lsplit<-j-i
        if(lsplit>(k/2)){
          lsplit<-k-lsplit
        }
        y<-paste(gsub(" ","",paste("[",s,","),perl = TRUE),gsub(" ","",paste("size=",lsplit,"]"),perl = TRUE),a[s])
        if(i==1){
          if(j==2){
            nsplit<-"1,"
          }
          if(j>2){
            nsplit<-c(1)
            if(j>3){
              for (ij in 2:(j-2)) {
                nsplit<-paste(nsplit,ij)
              }
            }
            nsplit<-paste(nsplit,gsub(" ","",paste((j-1),","),perl = TRUE))
          }
        }
        if(i!=1){
          nsplit<-c(1)
          if(i>2){
            for (ij in 2:(i-1)) {
              nsplit<-paste(nsplit,ij)
            }
          }
          if(j<k){
            for (ij in j:(k-1)) {
              nsplit<-paste(nsplit,ij)
            }
          }
          nsplit<-paste(nsplit,gsub(" ","",paste(k,","),perl = TRUE))
        }
        y<-paste(y,nsplit)
        x<-c(x,y)
      }
    }
  }
  x<-c(x,";")
  x<-c(x,"END; [Splits]")
  residue<-a_residue$deviance
  lsfit<-(1-residue/sum(p^2))*100
  return(list(x,lsfit))
}
#use ordering to write nexus file
draw_network_split_block<-function(a,A,taxaname=NULL){
  n<-length(a)
  z<-a
  for (i in 1:n) {
    x<-i
    while (a[i]!=i) {
      b<-A[a[i],]
      A[a[i],]<-A[x,]
      A[x,]<-b
      b<-A[,a[i]]
      A[,a[i]]<-A[,x]
      A[,x]<-b
      c<-a[a[i]]
      x<-a[a[i]]<-a[i]
      a[i]<-c
    }
  }
  np<-A[row(A)<col(A)]
  y<-make.nexus_split_block_order(z,np,n,taxaname)
  return(y)
}



#' Use an Input Tree for Lpnet
#'
#' Replacing the tree which constructed by the methods from \code{\link{lpnet}} to an arbitrary input tree which need to be input.
#'
#' @param M the distance matrix for construct network (the matrix should fit the triangle inequality and the diagonal should be 0).
#' @param tree the phylogenetic tree which class is \code{phylo} for lpnet algorithm (the \code{edge} block is necessary).
#' @param lp.package which package will used for Linear Programming, default is \code{Rglpk}, for a free R package;
#' \code{gurobi} for the gurobi package.
#' @param lp.type a character vector indicating the types of the objective variables. default is \code{B} for binary;
#' \code{I} for intrger; \code{C} for continuous; \code{NULL}, for ordinary.
#' @param  filename a character will be the naxus file's name, default is \code{lpnet.nex}.
#' @param  taxaname a character set of names for taxa, ordering is consist with original distance matrix \code{M}.
#'
#' @return The LSfit value.
#'
#' @importFrom ape nj
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom nnls nnls
#' @importFrom utils read.table
#' @importFrom utils write.table
#'
#' @examples
#' ### From Huson and Bryant (2006, Fig 4):
#' x <- c(14.06, 17.24, 20.5, 23.37, 17.43, 19.18, 18.48, 9.8, 13.06, 15.93, 15.65,
#'        17.4, 16.7, 6.74, 16.87, 16.59, 18.34, 17.64, 17.57, 17.29, 19.04, 18.34,
#'        17.6, 19.35, 21.21, 9.51, 11.37, 13.12)
#' M <- matrix(0, 8, 8)
#' M[row(M) > col(M)] <- x
#' M <- M + t(M)
#' njtree <- ape::nj(M)
#' taxaname <- c("A", "B", "C", "D", "E", "F", "G", "H")
#' lpnet.input.tree(M,
#'                  tree = njtree,
#'                  lp.package = "Rglpk",
#'                  lp.type = "B",
#'                  filename = "example.nex",
#'                  taxaname = taxaname)
#'
#' @export
lpnet.input.tree<-function(M,tree,lp.package="Rglpk",lp.type="B",filename="lpnet.nex",taxaname=NULL){
  n<-sqrt(length(M))#dim is n

  tr1<-tree

  tredge=tr1[["edge"]]#list of which taxa under a interior node
  taxa<-list()
  for (i in 1:n) {
    taxa[[i]]<-i
  }
  for (i in 1:(n-2)) {
    m<-2*n-1-i
    t1<-tredge[tredge[,1]==m][3]
    t2<-tredge[tredge[,1]==m][4]
    taxa[[m]]<-c(taxa[[t1]],taxa[[t2]])
  }

  #original ordering for the original tree
  original_ordering<-c()
  for (i in (1:(2*n-3))) {
    if(tredge[i,2]<=n){
      original_ordering<-c(original_ordering,tredge[i,2])
    }
  }

  if(lp.package=="Rglpk"){
    #calculate the parameter for lp
    lp_c<-make_lp_c(n,M,tredge,taxa)
    lp_matrix<-make_lp_matrix(n)
    lp_rhs<-make_lp_rhs(n)
    lp_dir<-make_lp_dir(n)

    if(is.null(lp.type)){#object variables are ordinary
      eg.rglpk<-Rglpk_solve_LP(obj = lp_c,
                               mat = lp_matrix,
                               rhs = lp_rhs,
                               dir = lp_dir,
                               max = TRUE)
    }else if(lp.type=="C"){#object variables are continuous
      eg.rglpk<-Rglpk_solve_LP(obj = lp_c,
                               mat = lp_matrix,
                               rhs = lp_rhs,
                               dir = lp_dir,
                               max = TRUE,
                               types = "C")
    }else if(lp.type=="B"){#object variables are binary
      eg.rglpk<-Rglpk_solve_LP(obj = lp_c,
                               mat = lp_matrix,
                               rhs = lp_rhs,
                               dir = lp_dir,
                               max = TRUE,
                               types = "B")
    }else if(lp.type=="I"){#object variables are integer
      eg.rglpk<-Rglpk_solve_LP(obj = lp_c,
                               mat = lp_matrix,
                               rhs = lp_rhs,
                               dir = lp_dir,
                               max = TRUE,
                               types = "I")
    }

    lpsolution<-eg.rglpk$solution
  }

  if(lp.package=="gurobi"){
    #calculate the parameter for lp
    model <- list()
    model$A <- make_lp_matrix(n)
    model$obj <- make_lp_c(n,M,tredge,taxa)
    model$modelsense <- 'max'
    model$rhs <- make_lp_rhs(n)
    model$sense <- make_lp_dir(n)
    params <- list(OutputFlag=0)
    if(is.null(lp.type)){#object variables are ordinary

    }else if(lp.type=="C"){#object variables are continuous
      model$vtype <- "C"
    }else if(lp.type=="B"){#object variables are binary
      model$vtype <- "B"
    }else if(lp.type=="I"){#object variables are integer
      model$vtype <- "I"
    }
    gurobi_lp<-gurobi(model,params)

    lpsolution<-gurobi_lp$x
  }

  lp_ordering<-lp_opt(original_ordering,lpsolution,tredge,taxa)#calculate lp net's circular ordering

  split_block_and_lsfit<-draw_network_split_block(lp_ordering,M,taxaname)

  split_block<-split_block_and_lsfit[[1]]

  write.table(split_block, file =filename, quote = FALSE, row.names = FALSE, col.names = FALSE)

  lsfit<-split_block_and_lsfit[[2]]

  return(lsfit)
}


