# Lpnet

## Package

A R package for constructing a circular split network from a distance matrix by Linear Programming. The functions in `lpnet` package are<br>
* `lpnet`<br>
* `read.nexus.distanceblock`<br>
* `read.nexus.taxablock`<br>

The main function is `lpnet`. `lpnet` constructs a planner network which has a circular ordering for a distance matrix, and writes a nexus file for [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)[1]. First construct a tree by different tree construct method (`neighbot-joining`[4], symmetry and not symmetry `neighbornet tree`, `UNJ`[3], `BioNJ`[2]) for the distance matrix.Then use Linear Programming(lp) to change the circular ordering. Then use Non-negative least squares(nnls) to calculate weights of splits which are consist with the lp net ordering. Finally, return a LSfit which is the least squares fit between the pairwise distances in the graph and the pairwise distances in the matrix. And write a nexus file with taxa block and splits block for [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) to see the circular network.<br>

The functions `read.nexus.distanceblock` and `read.nexus.taxablock` are tools to read the distance matrix and taxa names from a nex file instead of input them by hand. In addition the distance block in the nex file should be set as `diagonal` and `triangle=both`.

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
A character vector indicating the types of the objective variables. default is `NULL`, for ordinary; `B` for binary; `I` for intrger; `C` for continuous.<br>
* `filename`<br>
A character will be the naxus file's name, default is `lpnet.nex`.<br>
* `taxaname`<br>
A character set of names for taxa, ordering is consist with original distance matrix `M`.<br>
* `sequencelength`<br>
The sequence length of the data only used for BioNJ (default is 1 for only distance matrix).<br>

## Reference

[1][Huson, D. H., & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular biology and evolution, 23(2), 254-267.](https://academic.oup.com/mbe/article/23/2/254/1118872)<br>
[2][Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular biology and evolution, 14(7), 685-695.](https://academic.oup.com/mbe/article-abstract/14/7/685/1119804)<br>
[3][Gascuel, O. (1997). Concerning the NJ algorithm and its unweighted version, UNJ. Mathematical hierarchies and biology, 37, 149-171.](https://books.google.com/books?hl=zh-CN&lr=&id=stL67JmcWSkC&oi=fnd&pg=PA149&dq=Concerning+the+NJ+algorithm+and+its+unweighted+version,+UNJ&ots=WVM_Ligot1&sig=QaVyXPWnIs6R2090OTsmO41duBQ)<br>
[4][Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular biology and evolution, 4(4), 406-425.](https://academic.oup.com/mbe/article-abstract/4/4/406/1029664)
