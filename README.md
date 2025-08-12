# Branch Decomposition-Independent Edit Distances for Merge Trees

## Overview

This repository provides a Python library implementing the methods described in the following research papers:

"Branch Decomposition-Independent Edit Distances for Merge Trees."  
Florian Wetzels, Heike Leitte, and Christoph Garth.  
Computer Graphics Forum, 2022.  
DOI: [10.1111/cgf.14547](https://doi.org/https://doi.org/10.1111/cgf.14547)

"A Deformation-Based Edit Distance for Merge Trees"  
Florian Wetzels, Christoph Garth.  
TopoInVis 2022.  
DOI: [10.1109/TopoInVis57755.2022.00010](https://doi.org/10.1109/TopoInVis57755.2022.00010)

"Taming Horizontal Instability in Merge Trees: On the Computation of a Comprehensive Deformation-Based Edit Distance"
Florian Wetzels, Markus Anders, Christoph Garth.  
TopoInVis 2022.  
DOI: [0.1109/TopoInVis60193.2023.00015](https://doi.org/10.1109/TopoInVis60193.2023.00015)

#### Contents
The repository includes the library itself, as well as example scripts demonstrating how to use the library and replicate results from the papers.
For comparison with a method from the literature, the constrained edit distance for merge trees (see https://doi.org/10.1007/BF01975866) is also provided.
A more efficient C++ implementation of the path and branch mapping distance is available in a separate [repository]() or in the [TTK library]().

#### Original Source Code
For an archive of the original source code submitted with the first two papers on the branch and path mapping distance, see the [corresponding branch]().
This branch also contains an alternative C++ implementation, on which the benchmarks in the papers are based.
Furthermore, the original repositories contained datasets to showcase the distances.
These have been moved to separate data publications, but can also be found in the branch mentioned above.
The datasets are published here:

#### Datasets
"Vertical Instability Example — Four Clusters"  
Florian Wetzels, Heike Leitte, Christoph Garth.  
Zenodo, 2025.  
DOI: [10.5281/zenodo.16756130](https://doi.org/10.5281/zenodo.16756130)

"Vertical Instability Example — Outlier"  
Florian Wetzels, Heike Leitte, Christoph Garth.  
Zenodo, 2025.  
DOI: [10.5281/zenodo.16755706](https://doi.org/10.5281/zenodo.16755706)

"Horizontal Instability Example"  
Florian Wetzels, Markus Anders, Christoph Garth.  
Zenodo, 2025.  
DOI: [10.5281/zenodo.16758783](https://doi.org/10.5281/zenodo.16758783)

## Using the library

This library contains python implementations for four different tree edit distances:
- the classic constrained edit distance (details [here](https://doi.org/10.1007/BF01975866)) in `mted/constrained_edist_mt.py`;
- the branch mapping distance (details [here](https://doi.org/https://doi.org/10.1111/cgf.14547)) in `mted/branch_mapping_dist.py`;
- the path mapping distance (details [here](https://doi.org/10.1109/TopoInVis57755.2022.00010)) in `mted/path_mapping_dist.py`;
- the deformation-based edit distance (details [here](https://doi.org/10.1109/TopoInVis60193.2023.00015)) in `mted/deform_edist`.

For the constrained edit distance and the branch mapping distance, several base metrics to compare two branches are provided in `mted/baseMetrics.py`.
Usage is illustrated in the following examples:

```
from mted.baseMetrics import *
from mted.branch_mapping_dist import *

...

tree1_rooted = nx.dfs_tree(tree1, root1)
tree1_rooted.add_nodes_from((i, tree1.nodes[i]) for i in tree1_rooted.nodes)
tree2_rooted = nx.dfs_tree(tree2, root2)
tree2_rooted.add_nodes_from((i, tree2.nodes[i]) for i in tree2_rooted.nodes)

dist = branchMappingDistance(tree1_rooted,root1,tree2_rooted,root2,cost_wasserstein_branch,sqrt=True,traceback=False)
dist,match = branchMappingDistance(tree1_rooted,root1,tree2_rooted,root2,cost_wasserstein_branch,sqrt=True,traceback=True)
```

Here, `tree1` and `tree2` are the merge trees as networkx graphs with scalar values per node as well as birth and death values. The corresponding rooted trees `tree1_rooted` and `tree2_rooted` are computed using he DFS traversal provided by networkx. The indices of the merge tree roots are given as `root1` and `root2`. The cost function for branches `cost_wasserstein_branch` can be replaced by an arbitrary distance function that can be defined like in `mted/baseMetrics.py`. The return value is either just the distance or a tuple containing the distance and its corresponding matching, depending on whether the `traceback` flag is set. The `sqrt` flag determined whether the sum of all costs is returned or the square root of this sum.

The constrained merge tree edit distance can be used in exactly the same way:

```
from mted.baseMetrics import *
from mted.constrained_edist_mt import *

...

tree1_rooted = nx.dfs_tree(tree1, root1)
tree1_rooted.add_nodes_from((i, tree1.nodes[i]) for i in tree1_rooted.nodes)
tree2_rooted = nx.dfs_tree(tree2, root2)
tree2_rooted.add_nodes_from((i, tree2.nodes[i]) for i in tree2_rooted.nodes)

dist = editDistance_constrained(tree1_rooted,root1,tree2_rooted,root2,cost_wasserstein_branch_squared,sqrt=True,traceback=False)
dist,match = editDistance_constrained(tree1_rooted,root1,tree2_rooted,root2,cost_wasserstein_branch_squared,sqrt=True,traceback=True)
```

For the path mapping distance, usage is very similar, with some minor adaptations. The `sqrt` flag is not available, and the base metric is fixed as well:

```
from mted.baseMetrics import *
from mted.path_mapping_dist import *

...

tree1_rooted = nx.dfs_tree(tree1, root1)
tree1_rooted.add_nodes_from((i, tree1.nodes[i]) for i in tree1_rooted.nodes)
tree2_rooted = nx.dfs_tree(tree2, root2)
tree2_rooted.add_nodes_from((i, tree2.nodes[i]) for i in tree2_rooted.nodes)

dist = pathMappingDistance(tree1,root1,tree2,root2,traceback=False)
dist,match = pathMappingDistance(tree1,root1,tree1,root2,traceback=True)
```

For the deformation-based edit distance, the input trees are given as the original unrooted graphs, and no `traceback` flag is available. Otherwise, the interface is the same: 

```
from mted.baseMetrics import *
from mted.deform_edist import *

...

dist = rooted_unordered_deformation_distance(tree1,root1,tree2,root2)
```

For computation of the merge trees with TTK and conversion to networkx data structures, helper functions are provided in `mergetree/mergetree.py`:

```
from mergetree.mergetree import *

tree = getMergeTree(path,simplificationThreshold,fieldName,treeType)
```

Here, `path` is the path to a `.vti` file, `simplificationThreshold` it the relative threshold value (between 0 and 1) passed to TTK for topological simplification, `fieldname` is the name of the VTK array containing the scalar field, and `treeType` defines the type of merge tree (0 for join trees, 1 for split trees).

## Running the examples

Three example scripts are provided in the `examples/` directory. They contain code for replicating results from the papers.
In particular, distance matrices on three different synthetic benchmark datasets are computed to showcase the stability of the distances. The first two scripts reproduce results from the papers on the path mapping and branch mapping distance, only using the distances considered there. The last example reproduces results from the paper on the unconstrained deformation-based edit distance, including the MIP implementation of this distance. For this last example, the computation is parallelized (8 threads), due to the high computational load of the MIP implementation. The first two examples work sequentially on one thread.

The dataset for the examples can be retrieved from the data publications referenced above or corresponding git repositories.
In `scripts/get_data.sh`, a bash script retrieving the datasets through git is provided.

Based on these datasets, the example scripts work as follows:

- By running `python3 examples/distanceMatrixComparison_cluster.py`, we get three distance matrices on the dataset "Vertical Instability Example — Four Clusters": one is computed with the path mapping distance, one with the branch mapping distance, and one with the constrained edit distance. They are rendered as a clustermap using SciPy and seaborn. The results are stored in `clustermap_<distance>_cluster.svg`.

- By running `python3 examples/distanceMatrixComparison_outlier.py`, we get three distance matrices on the dataset "Vertical Instability Example — Outlier": one is computed with the path mapping distance, one with the branch mapping distance, and one with the constrained edit distance. They are rendered as a clustermap using SciPy and seaborn. The results are stored in `clustermap_<distance>_outlier.svg`.

- By running `python3 examples/distanceMatrixComparison_horizontal.py`, we get four distance matrices on the dataset "Horizontal Instability Example": one is computed with the path mapping distance, one with the branch mapping distance, one with the constrained edit distance, and one with the unconstrained deformation-based edit distance. They are rendered as a clustermap using SciPy and seaborn. The results are stored in `clustermap_<distance>_hs.svg`.

## Dependencies

The code in this repository depends on the following libraries that should be installed to run the examples:
- [TTK](https://topology-tool-kit.github.io/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [seaborn](https://seaborn.pydata.org/)

For the MIP implementation of the deformation-based edit distance, additional optimization libraries are required, including a valid license:
- [PuLP](https://coin-or.github.io/pulp/)
- [Gurobi](https://docs.gurobi.com/current/)

From a vanilla Ubuntu 22.04, the following commands install all the dependencies except the Gurobi license:

```
# install apt packages
apt install python3-sklearn
apt install wget

# install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

# install conda packages
conda create --name bdied
conda activate bdied
conda install conda-forge:vtk
conda install seaborn
conda install -c conda-forge topologytoolkit
conda install networkx
conda install -c gurobi gurobi
connda install scipy
conda install conda-forge::pulp
```
