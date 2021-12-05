# Lpnet

## Package

A R package which main part is constructing a circular split network from a distance matrix by Linear Programming. We also provide an alternative heuristic method which has lower accurancy but bigger taxa number limit. The functions in `lpnet` package are<br>
* `lpnet`<br>
* `lpnet.input.tree`<br>
* `read.nexus.distanceblock`<br>
* `read.nexus.taxablock`<br>
* `nnls.only.use.b`<br>
* `heuristic.method`<br>

The main function is `lpnet`. `lpnet` constructs a planar network which has a circular ordering for a distance matrix, and writes a nexus file for [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)[1]. First construct a tree by different tree construct method (`neighbot-joining`[4], symmetry and not symmetry `neighbornet tree`, `UNJ`[3], `BioNJ`[2]) for the distance matrix.Then use Linear Programming(lp) to change the circular ordering. Then use Non-negative least squares(nnls) to calculate weights of splits which are consist with the lp net ordering. Finally, return a LSfit which is the least squares fit between the pairwise distances in the graph and the pairwise distances in the matrix. And write a nexus file with taxa block and splits block for [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) to see the circular network.

The function `lpnet.input.tree` allows user to use an arbitrary phylogenetic binary tree to repalce the tree which constructed in `lpnet`. In addition the input tree's class should be `phylo` and the `edge` block is necessary.

The functions `read.nexus.distanceblock` and `read.nexus.taxablock` are tools to read the distance matrix and taxa names from a nex file instead of input them by hand. In addition the distance block in the nex file should be set as `diagonal` and `triangle=both`.

The function `nnls.only.use.b` is an intermediate step of `heuristic.method` and we provide this function for test. The nnls(Non-negative least squares) algorithm solves <a href="https://www.codecogs.com/eqnedit.php?latex=\min{\parallel&space;Ax-b\parallel_2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min{\parallel&space;Ax-b\parallel_2}" title="\min{\parallel Ax-b\parallel_2}" /></a> with the constraint <a href="https://www.codecogs.com/eqnedit.php?latex=x\geqslant&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x\geqslant&space;0" title="x\geqslant 0" /></a>. In this function, we construct the corresponding matrix `A` in `Fortran` instead of in `R` which will use less memory and run time. So we only need the distance matrix `M` which will be used to calculate the vector `b` directly and a circular ordering as the inputs for the `nnls.only.use.b` function.

The function `heuristic.method` is an alternative to `lpnet`. For the precess, the `heuristic.method` is similar with `lpnet`. The only different is the step which calculate the circular ordering from a tree. The heuristic method flip the edge from top to bottom which has the biggest `W` value(The `W` value is the average of all changed quartets weights over the changed quartets number if flip the edge). Finally, we cycle through an improvement step that checks all edges and flips which edges can improve the sum of all quartets until no improvement or the number of cycles reaches the loop limit. Actually, `lpnet` has higher accuracy than `heuristic.method`. So, We want to use the `heuristic.method` for more taxa number that `lpnat` can't handle yet.

## Installation

First, install the R package `devtools`.<br>
Then, install the `lpnet` package from github.<br>

    devtools::install_github('yukimayuli-gmz/lpnet',ref = "main")

## Calculate distance matrix from sequence alignment

The input of function `lpnet` is a distance matrix. We can use the software [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) to read a sequence alignment file(`.nex`, `.fasta`) and save as a nexus file which has a distance matrix(the distance block should be set as `diagonal` and `triangle=both`). In addition, [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) has a variety of methods and parameters for calculating distance matrix from sequence alignment. Then we can use function `read.nexus.distanceblock` to read the distance matrix in R.

## Parameters of `lpnet`

When use the `lpnet` function in R to construct a network, there some parameters can be set.<br>
* `M`<br>
The distance matrix for construct tree and network (the matrix should fit the triangle inequality and the diagonal should be 0).<br>
* `tree.method`<br>
Method for construct the original tree for lp, default is `unj`, for unweighted ntighbor joining tree; `nj` for neighbor joining tree; `nnet` for symmetry nnet tree; `nnetns` for no symmetry nnet tree; `BioNJ` for BioNJ tree.<br>
* `lp.package`<br>
Which package will used for Linear Programming, default is `Rglpk`, for a free R package; `gurobi` for the gurobi package.<br>
* `lp.type`<br>
A character vector indicating the types of the objective variables. default is `B`, for binary; `I` for intrger; `C` for continuous; `NULL` for ordinary.<br>
* `filename`<br>
A character will be the naxus file's name, default is `lpnet.nex`.<br>
* `taxaname`<br>
A character set of names for taxa, ordering is consist with original distance matrix `M`.<br>

## The LSfit value

The `lpnet` function will return the `LSfit` value in R. The `LSfit` will measure the fit between the pairwise distances in the constructed network and the pairwise distances in the origin matrix. And the `LSfit` is given by:<br>
<a href="https://www.codecogs.com/eqnedit.php?latex=(1-\frac{\sum_{ij}(d_{ij}-p_{ij})^{2}}{\sum_{ij}d_{ij}^{2}})*100." target="_blank"><img src="https://latex.codecogs.com/gif.latex?(1-\frac{\sum_{ij}(d_{ij}-p_{ij})^{2}}{\sum_{ij}d_{ij}^{2}})*100." title="(1-\frac{\sum_{ij}(d_{ij}-p_{ij})^{2}}{\sum_{ij}d_{ij}^{2}})*100." /></a>
Where the <a href="https://www.codecogs.com/eqnedit.php?latex=d_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?d_{ij}" title="d_{ij}" /></a> are pairwise distances in origin matrix and <a href="https://www.codecogs.com/eqnedit.php?latex=p_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_{ij}" title="p_{ij}" /></a> are pairwise distances in constructed network.<br>
For `NeighborNet` network's `LSfit`, we can make [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) show the `LSfit` in the `Preferences Window` by using `Edit`→`Preferences`→`Status Line` to check the `Show LSFit` box.

## Reference

[1][Huson, D. H., & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular biology and evolution, 23(2), 254-267.](https://academic.oup.com/mbe/article/23/2/254/1118872)<br>
[2][Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular biology and evolution, 14(7), 685-695.](https://academic.oup.com/mbe/article-abstract/14/7/685/1119804)<br>
[3][Gascuel, O. (1997). Concerning the NJ algorithm and its unweighted version, UNJ. Mathematical hierarchies and biology, 37, 149-171.](https://books.google.com/books?hl=zh-CN&lr=&id=stL67JmcWSkC&oi=fnd&pg=PA149&dq=Concerning+the+NJ+algorithm+and+its+unweighted+version,+UNJ&ots=WVM_Ligot1&sig=QaVyXPWnIs6R2090OTsmO41duBQ)<br>
[4][Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular biology and evolution, 4(4), 406-425.](https://academic.oup.com/mbe/article-abstract/4/4/406/1029664)
