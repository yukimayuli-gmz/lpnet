#nnet tree
#calculate Q criterion
calculate_Q<-function(M){
  d<-sqrt(length(M))
  Q<-matrix(0,d,d)
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      Q[i,j]<-Q[j,i]<-(d-2)*M[i,j]-sum(M[i,1:d])-sum(M[j,1:d])
    }
  }
  return(Q)
}
#calculate distace between u and other taxa (when select 1-1 or 2-2)
calculate_reduced_M_11_or_22_nnet<-function(M,i,j){
  d<-sqrt(length(M))
  N<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        N[a,b]<-N[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        N[(a-1),b]<-N[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        N[(a-2),b]<-N[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        N[a,(b-1)]<-N[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        N[(a-1),(b-1)]<-N[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        N[(a-2),(b-1)]<-N[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        N[a,(b-2)]<-N[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        N[(a-1),(b-2)]<-N[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        N[(a-2),(b-2)]<-N[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      N[c,(d-1)]<-N[(d-1),c]<-0.5*(M[c,i]+M[c,j])
    }
    if((c>i)&(c<j)){
      N[(c-1),(d-1)]<-N[(d-1),(c-1)]<-0.5*(M[c,i]+M[c,j])
    }
    if(c>j){
      N[(c-2),(d-1)]<-N[(d-1),(c-2)]<-0.5*(M[c,i]+M[c,j])
    }
  }
  return(N)
}
#calculate distace between u and other taxa (when select 1-2)
calculate_reduced_M_12_nnet<-function(M,i,j){
  d<-sqrt(length(M))
  N<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        N[a,b]<-N[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        N[(a-1),b]<-N[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        N[(a-2),b]<-N[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        N[a,(b-1)]<-N[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        N[(a-1),(b-1)]<-N[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        N[(a-2),(b-1)]<-N[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        N[a,(b-2)]<-N[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        N[(a-1),(b-2)]<-N[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        N[(a-2),(b-2)]<-N[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      N[c,(d-1)]<-N[(d-1),c]<-1/3*M[c,i]+2/3*M[c,j]
    }
    if((c>i)&(c<j)){
      N[(c-1),(d-1)]<-N[(d-1),(c-1)]<-1/3*M[c,i]+2/3*M[c,j]
    }
    if(c>j){
      N[(c-2),(d-1)]<-N[(d-1),(c-2)]<-1/3*M[c,i]+2/3*M[c,j]
    }
  }
  return(N)
}
#calculate distace between u and two reduced taxa
calculate_lengh<-function(M,i,j){
  d<-sqrt(length(M))
  t<-c()
  t[1]<-1
  t[2]<-1
  return(t)
}
#calculate label after reduce
reduce_label<-function(L,i,j,n){
  l<-length(L)
  L<-L[-i]
  L<-L[-(j-1)]
  L<-c(L,n+l-2)
  return(L)
}
#construct matrix of edge and label
NNetT<-function(M){
  n<-sqrt(length(M))
  tree<-list()
  T<-matrix(0,2*n-3,3)
  L<-c(1:n)
  for (m in 1:(n-2)) {
    Q<-calculate_Q(M)
    s<-n-m+1
    dia<-max(Q)
    for (i in 1:s) {
      Q[i,i]<-dia+1
    }
    t<-which(Q==min(Q))
    r<-c()
    r[1]<-t[1]%/%s
    r[2]<-t[1]%%s
    if(r[2]!=0){
      r[1]<-r[1]+1
    }
    if(r[2]==0){
      r[2]<-s
    }
    i<-min(r)
    j<-max(r)
    edge_ij<-calculate_lengh(M,i,j)
    T[2*m-1,1]<-2*n-m-1
    T[2*m-1,2]<-L[i]
    T[2*m-1,3]<-edge_ij[1]
    T[2*m,1]<-2*n-m-1
    T[2*m,2]<-L[j]
    T[2*m,3]<-edge_ij[2]
    if(((L[i]>n)&(L[j]>n))|((L[i]<=n)&(L[j]<=n))){
      M<-calculate_reduced_M_11_or_22_nnet(M,i,j)
    }
    if((L[i]<=n)&(L[j]>n)){
      M<-calculate_reduced_M_12_nnet(M,i,j)
    }
    L<-reduce_label(L,i,j,n)
  }
  T[2*n-3,1]<-L[2]
  T[2*n-3,2]<-L[1]
  T[2*n-3,3]<-1
  tree$edge<-cbind(T[,1],T[,2])
  tree$edge.length<-T[,3]
  tree$tip.label<-as.character(c(1:n))
  tree$Nnode<-n-2
  class(tree)<-"phylo"
  return(tree)
}

#nnet no symmetry tree
#calculate distance between u and other taxa (when select 1-1)
calculate_reduced_M_11_nnetns<-function(M,i,j){
  d<-sqrt(length(M))
  A<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        A[a,b]<-A[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        A[(a-1),b]<-A[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        A[(a-2),b]<-A[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        A[a,(b-1)]<-A[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        A[(a-1),(b-1)]<-A[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        A[(a-2),(b-1)]<-A[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        A[a,(b-2)]<-A[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        A[(a-1),(b-2)]<-A[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        A[(a-2),(b-2)]<-A[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      A[c,(d-1)]<-A[(d-1),c]<-0.5*(M[c,i]+M[c,j])
    }
    if((c>i)&(c<j)){
      A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-0.5*(M[c,i]+M[c,j])
    }
    if(c>j){
      A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-0.5*(M[c,i]+M[c,j])
    }
  }
  return(A)
}
#calculate distance between u and other taxa (when select 1-2)
calculate_reduced_M_12_nnetns<-function(M,i,j){
  d<-sqrt(length(M))
  A<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        A[a,b]<-A[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        A[(a-1),b]<-A[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        A[(a-2),b]<-A[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        A[a,(b-1)]<-A[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        A[(a-1),(b-1)]<-A[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        A[(a-2),(b-1)]<-A[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        A[a,(b-2)]<-A[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        A[(a-1),(b-2)]<-A[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        A[(a-2),(b-2)]<-A[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      A[c,(d-1)]<-A[(d-1),c]<-1/3*M[c,i]+2/3*M[c,j]
    }
    if((c>i)&(c<j)){
      A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-1/3*M[c,i]+2/3*M[c,j]
    }
    if(c>j){
      A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-1/3*M[c,i]+2/3*M[c,j]
    }
  }
  return(A)
}
#calculate distance between u and other taxa (when select 2-2) (j[3] is outside)
calculate_reduced_M_22_nnetns<-function(M,N,L1,L2,i,j,n){
  d<-sqrt(length(M))
  A<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        A[a,b]<-A[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        A[(a-1),b]<-A[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        A[(a-2),b]<-A[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        A[a,(b-1)]<-A[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        A[(a-1),(b-1)]<-A[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        A[(a-2),(b-1)]<-A[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        A[a,(b-2)]<-A[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        A[(a-1),(b-2)]<-A[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        A[(a-2),(b-2)]<-A[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  a0<-sample(c(1,2),1)
  if(a0==1){
    for (c in 1:d) {
      if(L1[c,1]<=n){
        if(c<i){
          A[c,(d-1)]<-A[(d-1),c]<-2/9*N[which(L2==L1[c,1]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,1]),which(L2==L1[j,3])]
        }
        if((c>i)&(c<j)){
          A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-2/9*N[which(L2==L1[c,1]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,1]),which(L2==L1[j,3])]
        }
        if(c>j){
          A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-2/9*N[which(L2==L1[c,1]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,1]),which(L2==L1[j,3])]
        }
      }
      if(L1[c,1]>n){
        if(c<i){
          A[c,(d-1)]<-A[(d-1),c]<-1/2*(2/9*N[which(L2==L1[c,2]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,2]),which(L2==L1[j,3])])+1/2*(2/9*N[which(L2==L1[c,3]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,3]),which(L2==L1[j,3])])
        }
        if((c>i)&(c<j)){
          A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-1/2*(2/9*N[which(L2==L1[c,2]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,2]),which(L2==L1[j,3])])+1/2*(2/9*N[which(L2==L1[c,3]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,3]),which(L2==L1[j,3])])
        }
        if(c>j){
          A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-1/2*(2/9*N[which(L2==L1[c,2]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,2]),which(L2==L1[j,3])])+1/2*(2/9*N[which(L2==L1[c,3]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,2])]+3/9*N[which(L2==L1[c,3]),which(L2==L1[j,3])])
        }
      }
    }
  }
  if(a0==2){
    for (c in 1:d) {
      if(L1[c,1]<=n){
        if(c<i){
          A[c,(d-1)]<-A[(d-1),c]<-3/9*N[which(L2==L1[c,1]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,3])]
        }
        if((c>i)&(c<j)){
          A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-3/9*N[which(L2==L1[c,1]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,3])]
        }
        if(c>j){
          A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-3/9*N[which(L2==L1[c,1]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,1]),which(L2==L1[j,3])]
        }
      }
      if(L1[c,1]>n){
        if(c<i){
          A[c,(d-1)]<-A[(d-1),c]<-1/2*(3/9*N[which(L2==L1[c,2]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,3])])+1/2*(3/9*N[which(L2==L1[c,3]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,3])])
        }
        if((c>i)&(c<j)){
          A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-1/2*(3/9*N[which(L2==L1[c,2]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,3])])+1/2*(3/9*N[which(L2==L1[c,3]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,3])])
        }
        if(c>j){
          A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-1/2*(3/9*N[which(L2==L1[c,2]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,2]),which(L2==L1[j,3])])+1/2*(3/9*N[which(L2==L1[c,3]),which(L2==L1[i,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[i,3])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,2])]+2/9*N[which(L2==L1[c,3]),which(L2==L1[j,3])])
        }
      }
    }
  }
  return(list(A,a0))
}
#calculate distance between u, v and other taxa (when select 1-2) (j[3] is outside)
calculate_reduced_N_12_nnetns<-function(N,L1,L2,i,j){
  d<-sqrt(length(N))
  x1<-which(L2==L1[i,1])
  y1<-which(L2==L1[j,2])
  z1<-which(L2==L1[j,3])
  t<-sort(c(x1,y1,z1))
  x<-t[1]
  y<-t[2]
  z<-t[3]
  A<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<x)&(b<x)){
        A[a,b]<-A[b,a]<-N[a,b]
      }
      if((a>x)&(a<y)&(b<x)){
        A[(a-1),b]<-A[b,(a-1)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b<x)){
        A[(a-2),b]<-A[b,(a-2)]<-N[a,b]
      }
      if((a>z)&(b<x)){
        A[(a-3),b]<-A[b,(a-3)]<-N[a,b]
      }
      if((a<x)&(b>x)&(b<y)){
        A[a,(b-1)]<-A[(b-1),a]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>x)&(b<y)){
        A[(a-1),(b-1)]<-A[(b-1),(a-1)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>x)&(b<y)){
        A[(a-2),(b-1)]<-A[(b-1),(a-2)]<-N[a,b]
      }
      if((a>z)&(b>x)&(b<y)){
        A[(a-3),(b-1)]<-A[(b-1),(a-3)]<-N[a,b]
      }
      if((a<x)&(b>y)&(b<z)){
        A[a,(b-2)]<-A[(b-2),a]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>y)&(b<z)){
        A[(a-1),(b-2)]<-A[(b-2),(a-1)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>y)&(b<z)){
        A[(a-2),(b-2)]<-A[(b-2),(a-2)]<-N[a,b]
      }
      if((a>z)&(b>y)&(b<z)){
        A[(a-3),(b-2)]<-A[(b-2),(a-3)]<-N[a,b]
      }
      if((a<x)&(b>z)){
        A[a,(b-3)]<-A[(b-3),a]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>z)){
        A[(a-1),(b-3)]<-A[(b-3),(a-1)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>z)){
        A[(a-2),(b-3)]<-A[(b-3),(a-2)]<-N[a,b]
      }
      if((a>z)&(b>z)){
        A[(a-3),(b-3)]<-A[(b-3),(a-3)]<-N[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<x){
      A[c,(d-2)]<-A[(d-2),c]<-2/3*N[c,x1]+1/3*N[c,y1]
      A[c,(d-1)]<-A[(d-1),c]<-1/3*N[c,y1]+2/3*N[c,z1]
    }
    if((c>x)&(c<y)){
      A[(c-1),(d-2)]<-A[(d-2),(c-1)]<-2/3*N[c,x1]+1/3*N[c,y1]
      A[(c-1),(d-1)]<-A[(d-1),(c-1)]<-1/3*N[c,y1]+2/3*N[c,z1]
    }
    if((c>y)&(c<z)){
      A[(c-2),(d-2)]<-A[(d-2),(c-2)]<-2/3*N[c,x1]+1/3*N[c,y1]
      A[(c-2),(d-1)]<-A[(d-1),(c-2)]<-1/3*N[c,y1]+2/3*N[c,z1]
    }
    if(c>z){
      A[(c-3),(d-2)]<-A[(d-2),(c-3)]<-2/3*N[c,x1]+1/3*N[c,y1]
      A[(c-3),(d-1)]<-A[(d-1),(c-3)]<-1/3*N[c,y1]+2/3*N[c,z1]
    }
  }
  A[(d-1),(d-2)]<-A[(d-2),(d-1)]<-1/3*N[x1,y1]+1/3*N[y1,z1]+1/3*N[x1,z1]
  return(A)
}
#calculate distance between u, v and other taxa (when select 2-2) (j[3] is outside)
calculate_reduced_N_22_nnetns<-function(N,L1,L2,i,j,a0){
  d<-sqrt(length(N))
  w1<-which(L2==L1[i,2])
  x1<-which(L2==L1[i,3])
  y1<-which(L2==L1[j,2])
  z1<-which(L2==L1[j,3])
  t<-sort(c(w1,x1,y1,z1))
  w<-t[1]
  x<-t[2]
  y<-t[3]
  z<-t[4]
  A<-matrix(0,d-2,d-2)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<w)&(b<w)){
        A[a,b]<-A[b,a]<-N[a,b]
      }
      if((a>w)&(a<x)&(b<w)){
        A[(a-1),b]<-A[b,(a-1)]<-N[a,b]
      }
      if((a>x)&(a<y)&(b<w)){
        A[(a-2),b]<-A[b,(a-2)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b<w)){
        A[(a-3),b]<-A[b,(a-3)]<-N[a,b]
      }
      if((a>z)&(b<w)){
        A[(a-4),b]<-A[b,(a-4)]<-N[a,b]
      }
      if((a<w)&(b>w)&(b<x)){
        A[a,(b-1)]<-A[(b-1),a]<-N[a,b]
      }
      if((a>w)&(a<x)&(b>w)&(b<x)){
        A[(a-1),(b-1)]<-A[(b-1),(a-1)]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>w)&(b<x)){
        A[(a-2),(b-1)]<-A[(b-1),(a-2)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>w)&(b<x)){
        A[(a-3),(b-1)]<-A[(b-1),(a-3)]<-N[a,b]
      }
      if((a>z)&(b>w)&(b<x)){
        A[(a-4),(b-1)]<-A[(b-1),(a-4)]<-N[a,b]
      }
      if((a<w)&(b>x)&(b<y)){
        A[a,(b-2)]<-A[(b-2),a]<-N[a,b]
      }
      if((a>w)&(a<x)&(b>x)&(b<y)){
        A[(a-1),(b-2)]<-A[(b-2),(a-1)]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>x)&(b<y)){
        A[(a-2),(b-2)]<-A[(b-2),(a-2)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>x)&(b<y)){
        A[(a-3),(b-2)]<-A[(b-2),(a-3)]<-N[a,b]
      }
      if((a>z)&(b>x)&(b<y)){
        A[(a-4),(b-2)]<-A[(b-2),(a-4)]<-N[a,b]
      }
      if((a<w)&(b>y)&(b<z)){
        A[a,(b-3)]<-A[(b-3),a]<-N[a,b]
      }
      if((a>w)&(a<x)&(b>y)&(b<z)){
        A[(a-1),(b-3)]<-A[(b-3),(a-1)]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>y)&(b<z)){
        A[(a-2),(b-3)]<-A[(b-3),(a-2)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>y)&(b<z)){
        A[(a-3),(b-3)]<-A[(b-3),(a-3)]<-N[a,b]
      }
      if((a>z)&(b>y)&(b<z)){
        A[(a-4),(b-3)]<-A[(b-3),(a-4)]<-N[a,b]
      }
      if((a<w)&(b>z)){
        A[a,(b-4)]<-A[(b-4),a]<-N[a,b]
      }
      if((a>w)&(a<x)&(b>z)){
        A[(a-1),(b-4)]<-A[(b-4),(a-1)]<-N[a,b]
      }
      if((a>x)&(a<y)&(b>z)){
        A[(a-2),(b-4)]<-A[(b-4),(a-2)]<-N[a,b]
      }
      if((a>y)&(a<z)&(b>z)){
        A[(a-3),(b-4)]<-A[(b-4),(a-3)]<-N[a,b]
      }
      if((a>z)&(b>z)){
        A[(a-4),(b-4)]<-A[(b-4),(a-4)]<-N[a,b]
      }
    }
  }
  if(a0==1){
    for (c in 1:d) {
      if(c<w){
        A[c,(d-3)]<-A[(d-3),c]<-4/9*N[c,w1]+3/9*N[c,x1]+2/9*N[c,y1]
        A[c,(d-2)]<-A[(d-2),c]<-1/9*N[c,x1]+2/9*N[c,y1]+6/9*N[c,z1]
      }
      if((c>w)&(c<x)){
        A[(c-1),(d-3)]<-A[(d-3),(c-1)]<-4/9*N[c,w1]+3/9*N[c,x1]+2/9*N[c,y1]
        A[(c-1),(d-2)]<-A[(d-2),(c-1)]<-1/9*N[c,x1]+2/9*N[c,y1]+6/9*N[c,z1]
      }
      if((c>x)&(c<y)){
        A[(c-2),(d-3)]<-A[(d-3),(c-2)]<-4/9*N[c,w1]+3/9*N[c,x1]+2/9*N[c,y1]
        A[(c-2),(d-2)]<-A[(d-2),(c-2)]<-1/9*N[c,x1]+2/9*N[c,y1]+6/9*N[c,z1]
      }
      if((c>y)&(c<z)){
        A[(c-3),(d-3)]<-A[(d-3),(c-3)]<-4/9*N[c,w1]+3/9*N[c,x1]+2/9*N[c,y1]
        A[(c-3),(d-2)]<-A[(d-2),(c-3)]<-1/9*N[c,x1]+2/9*N[c,y1]+6/9*N[c,z1]
      }
      if(c>z){
        A[(c-4),(d-3)]<-A[(d-3),(c-4)]<-4/9*N[c,w1]+3/9*N[c,x1]+2/9*N[c,y1]
        A[(c-4),(d-2)]<-A[(d-2),(c-4)]<-1/9*N[c,x1]+2/9*N[c,y1]+6/9*N[c,z1]
      }
    }
    A[(d-2),(d-3)]<-A[(d-3),(d-2)]<-1/9*N[w1,x1]+1/9*N[w1,y1]+1/9*N[x1,y1]+2/9*N[w1,z1]+2/9*N[x1,z1]+2/9*N[y1,z1]
  }
  if(a0==2){
    for (c in 1:d) {
      if(c<w){
        A[c,(d-3)]<-A[(d-3),c]<-6/9*N[c,w1]+2/9*N[c,x1]+1/9*N[c,y1]
        A[c,(d-2)]<-A[(d-2),c]<-2/9*N[c,x1]+3/9*N[c,y1]+4/9*N[c,z1]
      }
      if((c>w)&(c<x)){
        A[(c-1),(d-3)]<-A[(d-3),(c-1)]<-6/9*N[c,w1]+2/9*N[c,x1]+1/9*N[c,y1]
        A[(c-1),(d-2)]<-A[(d-2),(c-1)]<-2/9*N[c,x1]+3/9*N[c,y1]+4/9*N[c,z1]
      }
      if((c>x)&(c<y)){
        A[(c-2),(d-3)]<-A[(d-3),(c-2)]<-6/9*N[c,w1]+2/9*N[c,x1]+1/9*N[c,y1]
        A[(c-2),(d-2)]<-A[(d-2),(c-2)]<-2/9*N[c,x1]+3/9*N[c,y1]+4/9*N[c,z1]
      }
      if((c>y)&(c<z)){
        A[(c-3),(d-3)]<-A[(d-3),(c-3)]<-6/9*N[c,w1]+2/9*N[c,x1]+1/9*N[c,y1]
        A[(c-3),(d-2)]<-A[(d-2),(c-3)]<-2/9*N[c,x1]+3/9*N[c,y1]+4/9*N[c,z1]
      }
      if(c>z){
        A[(c-4),(d-3)]<-A[(d-3),(c-4)]<-6/9*N[c,w1]+2/9*N[c,x1]+1/9*N[c,y1]
        A[(c-4),(d-2)]<-A[(d-2),(c-4)]<-2/9*N[c,x1]+3/9*N[c,y1]+4/9*N[c,z1]
      }
    }
    A[(d-2),(d-3)]<-A[(d-3),(d-2)]<-2/9*N[w1,x1]+2/9*N[w1,y1]+1/9*N[x1,y1]+2/9*N[w1,z1]+1/9*N[x1,z1]+1/9*N[y1,z1]
  }
  return(A)
}
#calculate label after reduce (label_L1(1-1))
reduce_label_L1_11<-function(L1,L2,i,j,n){
  l<-length(L1[,1])
  L1<-rbind(L1,c(n+l-2,L1[i,1],L1[j,1]))
  L1<-L1[-i,]
  L1<-L1[-(j-1),]
  return(L1)
}
#calculate label after reduce (label_L1(1-2))
reduce_label_L1_12<-function(L1,L2,i,j,n){
  l<-length(L1[,1])
  L1<-rbind(L1,c(n+l-2,2*n+l-2,3*n+l-2))
  L1<-L1[-i,]
  L1<-L1[-(j-1),]
  return(L1)
}
#calculate label after reduce (label_L1(2-2))
reduce_label_L1_22<-function(L1,L2,i,j,n){
  l<-length(L1[,1])
  L1<-rbind(L1,c(n+l-2,2*n+l-2,3*n+l-2))
  L1<-L1[-i,]
  L1<-L1[-(j-1),]
  return(L1)
}
#calculate label after reduce (label_L2(1-2))
reduce_label_L2_12<-function(L1,L2,i,j,n){
  l<-length(L1[,1])
  x1<-which(L2==L1[i,1])
  y1<-which(L2==L1[j,2])
  z1<-which(L2==L1[j,3])
  t<-sort(c(x1,y1,z1))
  x<-t[1]
  y<-t[2]
  z<-t[3]
  L2<-L2[-x]
  L2<-L2[-(y-1)]
  L2<-L2[-(z-2)]
  L2<-c(L2,2*n+l-2,3*n+l-2)
  return(L2)
}
#calculate label after reduce (label_L2(2-2))
reduce_label_L2_22<-function(L1,L2,i,j,n){
  l<-length(L1[,1])
  w1<-which(L2==L1[i,2])
  x1<-which(L2==L1[i,3])
  y1<-which(L2==L1[j,2])
  z1<-which(L2==L1[j,3])
  t<-sort(c(w1,x1,y1,z1))
  w<-t[1]
  x<-t[2]
  y<-t[3]
  z<-t[4]
  L2<-L2[-w]
  L2<-L2[-(x-1)]
  L2<-L2[-(y-2)]
  L2<-L2[-(z-3)]
  L2<-c(L2,2*n+l-2,3*n+l-2)
  return(L2)
}
#construct matrix of edge and label
NNetT_no_sym<-function(M){
  n<-sqrt(length(M))
  N<-M
  tree<-list()
  T<-matrix(0,2*n-3,2)
  L1<-matrix(0,nrow = n,ncol = 3)
  L1[,1]<-c(1:n)
  L2<-c(1:n)
  for (m in 1:(n-2)) {
    Q<-calculate_Q(M)
    s<-n-m+1
    dia<-max(Q)
    for (i in 1:s) {
      Q[i,i]<-dia+1
    }
    t<-which(Q==min(Q))
    r<-c()
    r[1]<-t[1]%/%s
    r[2]<-t[1]%%s
    if(r[2]!=0){
      r[1]<-r[1]+1
    }
    if(r[2]==0){
      r[2]<-s
    }
    i<-min(r)
    j<-max(r)
    T[2*m-1,1]<-2*n-m-1
    T[2*m-1,2]<-L1[i,1]
    T[2*m,1]<-2*n-m-1
    T[2*m,2]<-L1[j,1]
    a<-0
    if(a==0){
      if((L1[i,1]<=n)&(L1[j,1]<=n)){
        M<-calculate_reduced_M_11_nnetns(M,i,j)
        L1<-reduce_label_L1_11(L1,L2,i,j,n)
        a<-1
      }
    }
    if(a==0){
      if((L1[i,1]<=n)&(L1[j,1]>n)){
        P<-calculate_Q(N)
        x1<-which(L2==L1[i,1])
        y1<-which(L2==L1[j,2])
        z1<-which(L2==L1[j,3])
        if(P[x1,y1]>P[x1,z1]){
          p<-L1[j,3]
          L1[j,3]<-L1[j,2]
          L1[j,2]<-p
        }
        M<-calculate_reduced_M_12_nnetns(M,i,j)
        N<-calculate_reduced_N_12_nnetns(N,L1,L2,i,j)
        L2<-reduce_label_L2_12(L1,L2,i,j,n)
        L1<-reduce_label_L1_12(L1,L2,i,j,n)
        a<-1
      }
    }
    if(a==0){
      if((L1[i,1]>n)&(L1[j,1]>n)){
        P<-calculate_Q(N)
        w1<-which(L2==L1[i,2])
        x1<-which(L2==L1[i,3])
        y1<-which(L2==L1[j,2])
        z1<-which(L2==L1[j,3])
        q<-order(c(P[w1,y1],P[w1,z1],P[x1,y1],P[x1,z1]))
        if(q[1]==1){
          p<-L1[i,3]
          L1[i,3]<-L1[i,2]
          L1[i,2]<-p
        }
        if(q[1]==2){
          p<-L1[i,3]
          L1[i,3]<-L1[i,2]
          L1[i,2]<-p
          p<-L1[j,3]
          L1[j,3]<-L1[j,2]
          L1[j,2]<-p
        }
        if(q[1]==4){
          p<-L1[j,3]
          L1[j,3]<-L1[j,2]
          L1[j,2]<-p
        }
        M0<-calculate_reduced_M_22_nnetns(M,N,L1,L2,i,j,n)
        M<-M0[[1]]
        N<-calculate_reduced_N_22_nnetns(N,L1,L2,i,j,M0[[2]])
        L2<-reduce_label_L2_22(L1,L2,i,j,n)
        L1<-reduce_label_L1_22(L1,L2,i,j,n)
        a<-1
      }
    }
  }
  T[2*n-3,1]<-L1[2,1]
  T[2*n-3,2]<-L1[1,1]
  tree$edge<-cbind(T[,1],T[,2])
  tree$edge.length<-seq(1,1,length.out = 2*n-3)
  tree$tip.label<-as.character(c(1:n))
  tree$Nnode<-n-2
  class(tree)<-"phylo"
  return(tree)
}

#nnet balance tree
#calculate distance between u and other taxa
calculate_reduced_M_balance<-function(M,i,j,n1,n2){
  d<-sqrt(length(M))
  N<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        N[a,b]<-N[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        N[(a-1),b]<-N[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        N[(a-2),b]<-N[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        N[a,(b-1)]<-N[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        N[(a-1),(b-1)]<-N[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        N[(a-2),(b-1)]<-N[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        N[a,(b-2)]<-N[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        N[(a-1),(b-2)]<-N[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        N[(a-2),(b-2)]<-N[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      N[c,(d-1)]<-N[(d-1),c]<-n1/(n1+n2)*M[c,i]+n2/(n1+n2)*M[c,j]
    }
    if((c>i)&(c<j)){
      N[(c-1),(d-1)]<-N[(d-1),(c-1)]<-n1/(n1+n2)*M[c,i]+n2/(n1+n2)*M[c,j]
    }
    if(c>j){
      N[(c-2),(d-1)]<-N[(d-1),(c-2)]<-n1/(n1+n2)*M[c,i]+n2/(n1+n2)*M[c,j]
    }
  }
  return(N)
}
#construct matrix of edge and label
BalanceT<-function(M){
  n<-sqrt(length(M))
  tree<-list()
  T<-matrix(0,2*n-3,3)
  L<-c(1:n)
  for (m in 1:(n-2)) {
    Q<-calculate_Q(M)
    s<-n-m+1
    dia<-max(Q)
    for (i in 1:s) {
      Q[i,i]<-dia+1
    }
    t<-which(Q==min(Q))
    r<-c()
    r[1]<-t[1]%/%s
    r[2]<-t[1]%%s
    if(r[2]!=0){
      r[1]<-r[1]+1
    }
    if(r[2]==0){
      r[2]<-s
    }
    i<-min(r)
    j<-max(r)
    n1<-L[i]
    n2<-L[j]
    edge_ij<-calculate_lengh(M,i,j)
    T[2*m-1,1]<-2*n-m-1
    T[2*m-1,2]<-L[i]
    T[2*m-1,3]<-edge_ij[1]
    T[2*m,1]<-2*n-m-1
    T[2*m,2]<-L[j]
    T[2*m,3]<-edge_ij[2]
    M<-calculate_reduced_M_balance(M,i,j,S[n1],S[n2])
    L<-reduce_label(L,i,j,n)
    S[L[length(L)]]<-S[n1]+S[n2]
  }
  T[2*n-3,1]<-L[2]
  T[2*n-3,2]<-L[1]
  T[2*n-3,3]<-1
  tree$edge<-cbind(T[,1],T[,2])
  tree$edge.length<-T[,3]
  tree$tip.label<-as.character(c(1:n))
  tree$Nnode<-n-2
  class(tree)<-"phylo"
  return(tree)
}

#BioNJ tree
#calculate distace between u and other taxa
calculate_reduced_M_BioNJ<-function(M,i,j,lamda,edge_ij){
  d<-sqrt(length(M))
  N<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        N[a,b]<-N[b,a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        N[(a-1),b]<-N[b,(a-1)]<-M[a,b]
      }
      if((a>j)&(b<i)){
        N[(a-2),b]<-N[b,(a-2)]<-M[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        N[a,(b-1)]<-N[(b-1),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        N[(a-1),(b-1)]<-N[(b-1),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        N[(a-2),(b-1)]<-N[(b-1),(a-2)]<-M[a,b]
      }
      if((a<i)&(b>j)){
        N[a,(b-2)]<-N[(b-2),a]<-M[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        N[(a-1),(b-2)]<-N[(b-2),(a-1)]<-M[a,b]
      }
      if((a>j)&(b>j)){
        N[(a-2),(b-2)]<-N[(b-2),(a-2)]<-M[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      N[c,(d-1)]<-N[(d-1),c]<-lamda*(M[c,i]-edge_ij[1])+(1-lamda)*(M[c,j]-edge_ij[2])
    }
    if((c>i)&(c<j)){
      N[(c-1),(d-1)]<-N[(d-1),(c-1)]<-lamda*(M[c,i]-edge_ij[1])+(1-lamda)*(M[c,j]-edge_ij[2])
    }
    if(c>j){
      N[(c-2),(d-1)]<-N[(d-1),(c-2)]<-lamda*(M[c,i]-edge_ij[1])+(1-lamda)*(M[c,j]-edge_ij[2])
    }
  }
  return(N)
}
#calculate the v value between u and other taxa
calculate_reduced_V_BioNJ<-function(V,i,j,lamda){
  d<-sqrt(length(V))
  N<-matrix(0,d-1,d-1)
  for (a in (1:(d-1))) {
    for (b in ((a+1):d)) {
      if((a<i)&(b<i)){
        N[a,b]<-N[b,a]<-V[a,b]
      }
      if((a>i)&(a<j)&(b<i)){
        N[(a-1),b]<-N[b,(a-1)]<-V[a,b]
      }
      if((a>j)&(b<i)){
        N[(a-2),b]<-N[b,(a-2)]<-V[a,b]
      }
      if((a<i)&(b>i)&(b<j)){
        N[a,(b-1)]<-N[(b-1),a]<-V[a,b]
      }
      if((a>i)&(a<j)&(b>i)&(b<j)){
        N[(a-1),(b-1)]<-N[(b-1),(a-1)]<-V[a,b]
      }
      if((a>j)&(b>i)&(b<j)){
        N[(a-2),(b-1)]<-N[(b-1),(a-2)]<-V[a,b]
      }
      if((a<i)&(b>j)){
        N[a,(b-2)]<-N[(b-2),a]<-V[a,b]
      }
      if((a>i)&(a<j)&(b>j)){
        N[(a-1),(b-2)]<-N[(b-2),(a-1)]<-V[a,b]
      }
      if((a>j)&(b>j)){
        N[(a-2),(b-2)]<-N[(b-2),(a-2)]<-V[a,b]
      }
    }
  }
  for (c in 1:d) {
    if(c<i){
      N[c,(d-1)]<-N[(d-1),c]<-lamda*V[c,i]+(1-lamda)*V[c,j]-lamda*(1-lamda)*V[i,j]
    }
    if((c>i)&(c<j)){
      N[(c-1),(d-1)]<-N[(d-1),(c-1)]<-lamda*V[c,i]+(1-lamda)*V[c,j]-lamda*(1-lamda)*V[i,j]
    }
    if(c>j){
      N[(c-2),(d-1)]<-N[(d-1),(c-2)]<-lamda*V[c,i]+(1-lamda)*V[c,j]-lamda*(1-lamda)*V[i,j]
    }
  }
  return(N)
}
#calculate distace between u and two reduced taxa
calculate_lengh_BioNJ<-function(M,i,j){
  d<-sqrt(length(M))
  t<-c()
  t[1]<-((d-2)*M[i,j]+sum(M[i,1:d])-sum(M[j,1:d]))/2/(d-2)
  t[2]<-((d-2)*M[i,j]-sum(M[i,1:d])+sum(M[j,1:d]))/2/(d-2)
  return(t)
}
#calculate lamda
calculate_lamda_BioNJ<-function(V,i,j){
  d<-sqrt(length(V))
  lamda<-0.5+(sum(V[j,1:d])-sum(V[i,1:d]))/2/(d-2)/V[i,j]
  return(lamda)
}
#construct matrix of edge and label
BioNJ<-function(M,sequence_length = 1){
  n<-sqrt(length(M))
  V<-M/sequence_length
  tree<-list()
  T<-matrix(0,2*n-3,3)
  L<-c(1:n)
  for (m in 1:(n-2)) {
    Q<-calculate_Q(M)
    s<-n-m+1
    dia<-max(Q)
    for (i in 1:s) {
      Q[i,i]<-dia+1
    }
    t<-which(Q==min(Q))
    r<-c()
    r[1]<-t[1]%/%s
    r[2]<-t[1]%%s
    if(r[2]!=0){
      r[1]<-r[1]+1
    }
    if(r[2]==0){
      r[2]<-s
    }
    i<-min(r)
    j<-max(r)
    edge_ij<-calculate_lengh_BioNJ(M,i,j)
    lamda<-calculate_lamda_BioNJ(V,i,j)
    T[2*m-1,1]<-2*n-m-1
    T[2*m-1,2]<-L[i]
    T[2*m-1,3]<-edge_ij[1]
    T[2*m,1]<-2*n-m-1
    T[2*m,2]<-L[j]
    T[2*m,3]<-edge_ij[2]
    M<-calculate_reduced_M_BioNJ(M,i,j,lamda,edge_ij)
    V<-calculate_reduced_V_BioNJ(V,i,j,lamda)
    L<-reduce_label(L,i,j,n)
  }
  T[2*n-3,1]<-L[2]
  T[2*n-3,2]<-L[1]
  T[2*n-3,3]<-M[1,2]
  tree$edge<-cbind(T[,1],T[,2])
  tree$edge.length<-T[,3]
  tree$tip.label<-as.character(c(1:n))
  tree$Nnode<-n-2
  class(tree)<-"phylo"
  return(tree)
}

#adjust the ordering for edge in tree matrix
adjust_edge<-function(a,tredge,n){
  y<-a
  A<-matrix(0,nrow = 2*n-3,ncol = 2)
  t<-1
  ly<-length(y)
  while (ly>0) {
    a<-0
    if(a==0){
      if(y[1]<=n){
        p<-y[1]
        y<-y[-1]
        A[t,]<-tredge[tredge[,2]==p]
        t<-t+1
        a<-1
      }
    }
    if(a==0){
      if(y[1]>n){
        p<-y[1]
        y<-y[-1]
        y<-c(tredge[tredge[,1]==p][3],tredge[tredge[,1]==p][4],y)
        A[t,]<-tredge[tredge[,2]==p]
        t<-t+1
        a<-1
      }
    }
    ly<-length(y)
  }
  return(A)
}

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



#' Consruct a circular network based on a tree by Linear Programming
#'
#' Construct a planner network which has a circular ordering for a distance matrix, and write a nexus file for SplitsTree4.
#' First construct a tree for the distance matrix.Then use Linear Programming(lp) to change the circular ordering.
#' The ordering have the biggest sum of quartets for all taxa is the lp net ordering.
#' Then use Non-negative least squares(nnls) to calculate weights of splits which are consist with the lp net ordering.
#' Finally, return a LSfit which is the least squares fit between the pairwise distances in the graph and the pairwise distances in the matrix.
#' And write a nexus file with taxa block and splits block for SplitsTree4 to see the circular network.
#'
#' @param M the distance matrix for construct tree and network;
#' (the matrix should fit the triangle inequality and the diagonal should be 0).
#' @param tree.method method for construct the original tree for lp, default is \code{unj}, for unweighted ntighbor joining tree;
#' \code{nj} for neighbor joining tree; \code{nnet} for symmetry nnet tree; \code{nnetns} for no symmetry nnet tree;
#' \code{BioNJ} for BioNJ tree.
#' @param lp.package which package will used for Linear Programming, default is \code{Rglpk}, for a free R package;
#' \code{gurobi} for the gurobi package.
#' @param lp.type a character vector indicating the types of the objective variables. default is \code{NULL}, for ordinary;
#' \code{B} for binary; \code{I} for intrger; \code{C} for continuous.
#' @param  filename a character will be the naxus file's name, default is \code{lpnet.net}
#' @param  taxaname a character set of names for taxa, ordering is consist with original distance matrix \code{M}.
#' @param  sequencelength the sequence length of the data only used for BioNJ (default is 1 for only distance matrix).
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
#' x <- c(14.06, 17.24, 20.05, 23.37, 17.43, 19.18, 18.48, 9.8, 13.06, 15.93, 15.65,
#'        17.4, 16.7, 6.74, 16.87, 16.59, 18.34, 17.64, 17.57, 17.29, 19.04, 18.34,
#'        17.6, 19.35, 21.21, 9.51, 11.37, 13.12)
#' M <- matrix(0, 8, 8)
#' M[row(M) > col(M)] <- x
#' M <- M + t(M)
#' taxaname<-c("A", "B", "C", "D", "E", "F", "G", "H")
#' lpnet(M,
#'       tree.method = "nj",
#'       lp.package = "Rglpk",
#'       lp.type = NULL,
#'       filename = "example.nex",
#'       taxaname = taxaname)
#'
#' @export
lpnet<-function(M,tree.method="unj",lp.package="Rglpk",lp.type=NULL,filename="lpnet.nex",taxaname=NULL,sequencelength=1){
  n<-sqrt(length(M))#dim is n

  if(tree.method=="nj"){
    tr1 <- ape::nj(M)#tr1 is neighbor joining tree
  }
  if(tree.method=="nnet"){
    tr1 <- NNetT(M)#tr1 is nnet tree
    tredge=tr1[["edge"]]
    y<-c(tredge[tredge[,1]==n+1][4],tredge[tredge[,1]==n+1][5],tredge[tredge[,1]==n+1][6])
    tr1$edge<-adjust_edge(y,tredge,n)
  }
  if(tree.method=="nnetns"){
    tr1 <- NNetT_no_sym(M)#tr1 is nnet no symmetry tree
    tredge=tr1[["edge"]]
    y<-c(tredge[tredge[,1]==n+1][4],tredge[tredge[,1]==n+1][5],tredge[tredge[,1]==n+1][6])
    tr1$edge<-adjust_edge(y,tredge,n)
  }
  if(tree.method=="unj"){
    tr1 <- BalanceT(M)#tr1 is unweighted neighbor joining tree
    tredge=tr1[["edge"]]
    y<-c(tredge[tredge[,1]==n+1][4],tredge[tredge[,1]==n+1][5],tredge[tredge[,1]==n+1][6])
    tr1$edge<-adjust_edge(y,tredge,n)
  }
  if(tree.method=="BioNJ"){
    tr1 <- BioNJ(M,sequencelength)#tr1 is BioNJ tree
    tredge=tr1[["edge"]]
    y<-c(tredge[tredge[,1]==n+1][4],tredge[tredge[,1]==n+1][5],tredge[tredge[,1]==n+1][6])
    tr1$edge<-adjust_edge(y,tredge,n)
  }

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


