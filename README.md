# lpnet {lpnet}
=======================
## Consruct a circular network based on a tree by Linear Programming
----------------------------------------------------------------------
### Description
Construct a planner network which has a circular ordering for a distrance matrix, and write a nexus file for SplitsTree4. First construct a tree for the distance matrix.Then use Linear Programming(lp) to change the circular ordering. The ordering have the biggest sum of quartets for all taxa is the lp net ordering. Then use Non-negative least squares(nnls) to calculate weights of splits which are consist with the lp net ordering. Finally, write a nexus file with taxa block and splits block for SplitsTree4 to see the circular network.

### Usage<br>
lpnet(<br>
  M,<br>
  tree.method = "nj",<br>
  lp.package = "Rglpk",<br>
  lp.type = NULL,<br>
  filename = "lpnet.nex",<br>
  taxaname = NULL<br>
)

### Arguments
#### M	
the distance matrix for construct tree and network; (the matrix should fit the triangle inequality and the diagonal should be 0).

#### tree.method	
method for construct the original tree for lp, default is nj, for ntighbor joining; nnet for nnet tree; nnetns for no symmetry nnet tree; balance for balance nnet tree; qj for quartet joining tree.

#### lp.package	
which package will used for Linear Programming, default is Rglpk, for a free R package; gurobi for the gurobi package.

#### lp.type	
a character vector indicating the types of the objective variables. default is NULL, for ordinary; B for binary; I for intrger; C for continuous.

#### filename	
a character will be the naxus file's name, default is lpnet.net

#### taxaname	
a character set of names for taxa, ordering is consist with original distance matrix M.

#### Value
None (invisible ‘NULL’).

### Examples
#### From Huson and Bryant (2006, Fig 4):
x <- c(14.06, 17.24, 20.05, 23.37, 17.43, 19.18, 18.48, 9.8, 13.06, 15.93, 15.65,
       17.4, 16.7, 6.74, 16.87, 16.59, 18.34, 17.64, 17.57, 17.29, 19.04, 18.34,
       17.6, 19.35, 21.21, 9.51, 11.37, 13.12)
M <- matrix(0, 8, 8)
M[row(M) > col(M)] <- x
M <- M + t(M)
taxaname<-c("A", "B", "C", "D", "E", "F", "G", "H")
lpnet(M,
      tree.method = "nj",
      lp.package = "Rglpk",
      lp.type = NULL,
      filename = "example.nex",
      taxaname = taxaname)
