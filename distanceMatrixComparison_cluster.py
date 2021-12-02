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
    path = "./test_datasets/branchVSconstrained_cluster/test_MB/test_MB_"+str(i)+".vti"
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

#distanceMatrix_branch = np.zeros((20,20))
distanceMatrix_branch_sqrt = np.zeros((20,20))
#distanceMatrix_constrained = np.zeros((20,20))
distanceMatrix_constrained_sqrt = np.zeros((20,20))

for k in range(0,20):
    nodeScalars1,nodeLabels1,nodeTypes1,nodePositions1,upIDs1,downIDs1,regionSizes1,persistences1,children1,rootID1 = mergeTrees[k]
    nodeLabels1_zipped = list(zip(nodeTypes1,nodeLabels1))
    for l in range(0,20):
        nodeScalars2,nodeLabels2,nodeTypes2,nodePositions2,upIDs2,downIDs2,regionSizes2,persistences2,children2,rootID2 = mergeTrees[l]
        nodeLabels2_zipped = list(zip(nodeTypes2,nodeLabels2))
        
        print(k,l)
        
        #dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_wasserstein_branch,False)
        sqrtDist = math.sqrt(branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_wasserstein_branch_squared,False))
        #dist2 = editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_wasserstein_split,False)
        sqrtDist2 = math.sqrt(editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_wasserstein_split_squared,False))
        
#        dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_Linf_branch,False)
#        sqrtDist = math.sqrt(branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_Linf_branch_squared,False))
#        dist2 = editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_Linf_split,False)
#        sqrtDist2 = math.sqrt(editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_Linf_split_squared,False))

#        dist = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_overhang_branch,False)
#        sqrtDist = math.sqrt(branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_overhang_branch_squared,False))
#        dist2 = editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_overhang_split,False)
#        sqrtDist2 = math.sqrt(editDistance_constrained(nodeLabels1_zipped,children1,rootID1,nodeLabels2_zipped,children2,rootID2,cost_overhang_split_squared,False))
        
        print(k,l)
        print(k,l,sqrtDist)
        print(k,l,sqrtDist2)
        
        #distanceMatrix_branch[k,l] = dist
        #distanceMatrix_branch[l,k] = dist
        distanceMatrix_branch_sqrt[k,l] = sqrtDist
        #distanceMatrix_branch_sqrt[l,k] = sqrtDist
        #distanceMatrix_constrained[k,l] = dist2
        #distanceMatrix_constrained[l,k] = dist2
        distanceMatrix_constrained_sqrt[k,l] = sqrtDist2
        #distanceMatrix_constrained_sqrt[l,k] = sqrtDist2

#print(distanceMatrix)

dm_file = open("./test_datasets/branchVSconstrained_cluster/matrix_wasserstein_cpp.txt",'r')
dm_lines = dm_file.readlines()
dm_file.close()
dm_string = []
for line in dm_lines:
    dm_string.append(line.split(' '))

distanceMatrix_wasserstein = np.zeros((20,20))
for i in range(20):
    for j in range(20):
        distanceMatrix_wasserstein[i,j] = float(dm_string[i][j])
        
for i in range(20):
    for j in range(i,20):
        if(distanceMatrix_branch_sqrt[i,j]-distanceMatrix_branch_sqrt[j,i] != 0):
            distanceMatrix_branch_sqrt[i,j] = distanceMatrix_branch_sqrt[j,i]
        if(distanceMatrix_constrained_sqrt[i,j]-distanceMatrix_constrained_sqrt[j,i] != 0):
            distanceMatrix_constrained_sqrt[i,j] = distanceMatrix_constrained_sqrt[j,i]
        if(distanceMatrix_wasserstein[i,j]-distanceMatrix_wasserstein[j,i] != 0):
            print(distanceMatrix_wasserstein[i,j],distanceMatrix_wasserstein[j,i],distanceMatrix_wasserstein[i,j]-distanceMatrix_wasserstein[j,i])
            distanceMatrix_wasserstein[i,j] = distanceMatrix_wasserstein[j,i]

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_branch_sqrt), method='average')
cm = sns.clustermap(distanceMatrix_branch_sqrt, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_branch = cm.dendrogram_row.reordered_ind

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_constrained_sqrt), method='average')
cm = sns.clustermap(distanceMatrix_constrained_sqrt, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_constrained = cm.dendrogram_row.reordered_ind

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_wasserstein), method='average')
cm = sns.clustermap(distanceMatrix_wasserstein, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_wasserstein = cm.dendrogram_row.reordered_ind



fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_branch_sqrt[:, reorder_branch][reorder_branch], cmap='inferno_r', interpolation='nearest')#, vmin=0, vmax=maxV)
plt.axis('off')
#fig.colorbar(pos, ax=ax)
plt.colorbar(pos,ax=[ax],location='left',shrink=0.85)
fig.savefig("./clustermap_branch_cluster.svg",transparent=True,bbox_inches='tight',pad_inches=0)
#plt.show()

fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_constrained_sqrt[:, reorder_constrained][reorder_constrained], cmap='inferno_r', interpolation='nearest')#, vmin=0, vmax=maxV)
plt.axis('off')
#fig.colorbar(pos, ax=ax)
plt.colorbar(pos,ax=[ax],location='left',shrink=0.85)
fig.savefig("./clustermap_constrained_cluster.svg",transparent=True,bbox_inches='tight',pad_inches=0)
#plt.show()

fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_wasserstein[:, reorder_wasserstein][reorder_wasserstein], cmap='inferno_r', interpolation='nearest')#, vmin=0, vmax=maxV)
plt.axis('off')
#fig.colorbar(pos, ax=ax)
plt.colorbar(pos,ax=[ax],location='left',shrink=0.85)
fig.savefig("./clustermap_wasserstein_cluster.svg",transparent=True,bbox_inches='tight',pad_inches=0)


plt.show()

exit()
