# Branch Decomposition-Independent Edit Distances for Merge Trees

This repository contains a python library for computing branch decomposition-independent edit distances for merge trees as well as code examples and datasets to explain basic usage.

## Using the library

This library contains python implementations for two different tree edit distances, the classic constrained edit distance (see https://doi.org/10.1007/BF01975866) and the branch mapping distance in `mergeTreeEdit_branch.py` or `mergeTreeEdit_constrained.py`, as well as corresponding base metrics for merge trees in `baseMetrics.py`. Usage is illustrated in the following examples:

```
from baseMetrics import *
from mergeTreeEdit_branch import *

...

dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_wasserstein_branch,False)
dist,match = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_wasserstein_branch,True)
```

Here, `nodeScalars1[i]` and `nodeScalars2[i]` contain the scalar values of the node `i`, `children1[i]` and `children2[i]` contain the child nodes of node `i` and `rootID1`/`rootID2` are the indices of the two roots. `cost_wasserstein_branch` can be replaced by an arbitrary distance function that can be defined like in `baseMetrics.py`. The return value is either just the distance or a tuple containing the distance and its corresponding matching, depending on whether the traceback flag is set.

```
from baseMetrics import *
from mergeTreeEdit_constrained import *

...

dist = editDistance_constrained(nodeLabels1,children1,rootID1,nodeLabels2,children2,rootID2,cost_wasserstein_split_squared,False)
dist,match = editDistance_constrained(nodeLabels1,children1,rootID1,nodeLabels2,children2,rootID2,cost_wasserstein_split_squared,True)
```

Here, `nodeScalars1` and `nodeScalars2` contain arbitrary labels and everything else is defined as before.

## Running the examples

The four example files contain code for three different use cases: clustering, outlier detection and feature tracking.

The examples work as follows:

- By running `python3 distanceMatrixComparison_clustering.py`, we get three distance matrices for a synthetic example dataset (`test_datasets/branchVSconstrained_cluster/`): one is computed with the branch mapping distance, one with the constrained edit distance and one is a precomputed distance matrix using the TTKMergeTreeDistanceMatrix module (see https://doi.org/10.1109/TVCG.2021.3114839). They are rendered as a simple heatmap as well as a clustermap using SciPy.

- By running `python3 distanceMatrixComparison_outlier.py` or `python3 distanceMatrixComparison_outlier2.py`, we get three distance matrices for a synthetic example dataset (`test_datasets/branchVSconstrained_outlier/` or `test_datasets/branchVSconstrained_outlier2/`): one is computed with the branch mapping distance, one with the constrained edit distance and one is a precomputed distance matrix using the TTKMergeTreeDistanceMatrix module. They are rendered as a simple heatmap as well as a clustermap using SciPy.

- By running `python3 timeTrack_heatedCylinder2D_branch.py`, we get a vtk multiblock with feature tracking information attached on the heated cylinder dataset (`heatedCylinder2d.cdb/`). The vtk file is written to `output/`.

## Datasets and C++ version

The four datasets corresponding to the python examples can be found in the directories `heatedCylinder2d.cdb`, `test_datasets/branchVSconstrained_cluster`, `test_datasets/branchVSconstrained_outlier` and `test_datasets/branchVSconstrained_outlier2`. We also provide two noisy versions of the outlier ensemble: one in `test_datasets/branchVSconstrained_outlier2` and `test_datasets/branchVSconstrained_outlier_noisy`. All datasets are provided in VTK file formats. The two noisy datasets have no corresponding python examples, as computation times would be too high. The should only be used with the C++ version of the library (see below for details).

We also provide a more efficient C++ implementation of the two distance functions in the `cpp/` directory. The interface is similar to the python version. The C++ implementation does not yet output the matching, but only the distance, and does not allow for a base metric change. The Wasserstein base metric has been hardcoded for both distances.

## Dependencies

The code in this repository depends on the following libraries that should be installed to run the examples:
- [TTK](https://topology-tool-kit.github.io/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
