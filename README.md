# Lpnet

## Package

A R package for constructing a circular split network from a distance matrix by Linear Programming. The functions in `lpnet` package are<br>
* `lpnet`<br>
* `read.nexus.distanceblock`<br>
* `read.nexus.taxablock`<br>

The main function is `lpnet`. `lpnet` constructs a planner network which has a circular ordering for a distance matrix, and writes a nexus file for `SplitsTree4`. First construct a tree by different tree construct method (`neighbot-joining`, symmetry and not symmetry `neighbornet tree`, `UNJ`, `BioNJ`) for the distance matrix.Then use Linear Programming(lp) to change the circular ordering. Then use Non-negative least squares(nnls) to calculate weights of splits which are consist with the lp net ordering. Finally, return a LSfit which is the least squares fit between the pairwise distances in the graph and the pairwise distances in the matrix. And write a nexus file with taxa block and splits block for [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) to see the circular network.<br>

The functions `read.nexus.distanceblock` and `read.nexus.taxablock` are tools to read the distance matrix and taxa names from a nex file instead of input them by hand. In addition the distance block in the nex file should be set as `diagonal` and `triangle=both`.

## Installation

First, install the R package `devtools`.<br>
Then, install the `lpnet` package from github.<br>

    devtools::install_github('yukimayuli-gmz/lpnet',ref = "main")

## Calculate distance matrix from sequence alignment

The input of function `lpnet` is a distance matrix. We can use the software [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) to read a sequence alignment file(`.nex`, `.fasta`) and save as a nexus file which has a distance matrix(the distance block should be set as `diagonal` and `triangle=both`). In addition, [SplitsTree4](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) has a variety of methods and parameters for calculating distance matrix from sequence alignment. Then we can use function `read.nexus.distanceblock` to read the distance matrix in R.

## Reference

[Bryant, D., & Moulton, V. (2004). Neighbor-net: an agglomerative method for the construction of phylogenetic networks. Molecular biology and evolution, 21(2), 255-265.](https://academic.oup.com/mbe/article-abstract/21/2/255/1187993)<br>
[Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular biology and evolution, 14(7), 685-695.](https://academic.oup.com/mbe/article-abstract/14/7/685/1119804)<br>
[Gascuel, O. (1997). Concerning the NJ algorithm and its unweighted version, UNJ. Mathematical hierarchies and biology, 37, 149-171.](https://books.google.com/books?hl=zh-CN&lr=&id=stL67JmcWSkC&oi=fnd&pg=PA149&dq=Concerning+the+NJ+algorithm+and+its+unweighted+version,+UNJ&ots=WVM_Ligot1&sig=QaVyXPWnIs6R2090OTsmO41duBQ)<br>
[Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular biology and evolution, 4(4), 406-425.](https://academic.oup.com/mbe/article-abstract/4/4/406/1029664)
