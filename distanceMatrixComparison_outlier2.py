import math
import threading
import multiprocessing
import time

import vtk
import topologytoolkit as ttk
    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

from contourMergeTrees_helpers import *
from baseMetrics import *
from mergeTreeEdit_branch import *
from mergeTreeEdit_constrained import *

mergeTrees = []
for i in range(0,20):
    path = "./test_datasets/branchVSconstrained_outlier2/test_MB/test_MB_"+str(i)+".vti"
    simplificationThreshold = 0.0
    fieldName = "test"
    treeType = int(1)
    nodeScalars,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences = getMergeTree(path,simplificationThreshold,fieldName,treeType)
    
    rootID = -1
    children = []
    for i in range(len(nodeScalars)):
        if(treeType==0 and nodeTypes[i]==3):
            rootID = i
        if(treeType==1 and nodeTypes[i]==0):
            rootID = i
        children.append([])
    for i in range(len(upIDs)):
        downId = downIDs[i]
        upId = upIDs[i]
        children[upId].append(downId)
            
    bD,bL = getBranchDecomposition(children,nodeScalars,rootID)

    nodes = []
    persistences = dict()
    for i in range(len(bD)):
        nodes.append(bL[bD[i]])
        persistences[bD[i]] = bL[bD[i]][0] - bL[bD[i]][1]
    
    mergeTrees.append((nodeScalars,nodes,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences,children,rootID))
    
#print(mergeTrees)

distanceMatrix_branch = np.zeros((20,20))
distanceMatrix_constrained = np.zeros((20,20))

for k in range(0,20):
    nodeScalars1,nodeLabels1,nodeTypes1,nodePositions1,upIDs1,downIDs1,regionSizes1,persistences1,children1,rootID1 = mergeTrees[k]
    nodeLabels1_zipped = list(zip(nodeTypes1,nodeLabels1))
    for l in range(0,20):
        nodeScalars2,nodeLabels2,nodeTypes2,nodePositions2,upIDs2,downIDs2,regionSizes2,persistences2,children2,rootID2 = mergeTrees[l]
        nodeLabels2_zipped = list(zip(nodeTypes2,nodeLabels2))
        
        print(k,l)
        
        dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_wasserstein_branch,False)
        dist2 = editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_wasserstein_split,False)
        
#        dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_Linf_branch,False)
#        dist2 = editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_Linf_split,False)

#        dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_overhang_branch,False)
#        dist2 = editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_overhang_split,False)
        
        print(k,l)
        print(k,l,dist)
        print(k,l,dist2)
        
        distanceMatrix_branch[k,l] = dist
        #distanceMatrix_branch[l,k] = dist
        distanceMatrix_constrained[k,l] = dist2
        #distanceMatrix_constrained[l,k] = dist2

#print(distanceMatrix)

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_branch), method='average')
cm = sns.clustermap(distanceMatrix_branch, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_branch = cm.dendrogram_row.reordered_ind

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_constrained), method='average')
cm = sns.clustermap(distanceMatrix_constrained, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_constrained = cm.dendrogram_row.reordered_ind



fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_branch[:, reorder_branch][reorder_branch], cmap='inferno_r', interpolation='nearest')#, vmin=0, vmax=maxV)
plt.axis('off')
fig.colorbar(pos, ax=ax)
fig.savefig("./clustermap_branch_outlier2.svg",transparent=True,bbox_inches='tight',pad_inches=0)
#plt.show()

fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_constrained[:, reorder_constrained][reorder_constrained], cmap='inferno_r', interpolation='nearest')#, vmin=0, vmax=maxV)
plt.axis('off')
fig.colorbar(pos, ax=ax)
fig.savefig("./clustermap_constrained_outlier2.svg",transparent=True,bbox_inches='tight',pad_inches=0)
#plt.show()


plt.show()

exit()
