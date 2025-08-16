import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

from mted.baseMetrics import *
from mted.path_mapping_dist import *
from mted.branch_mapping_dist import *
from mted.constrained_edist_mt import *
from mergetree.mergetree import *

# compute merge trees for each scalar field

mergeTrees = []
for i in range(0,20):
    path = "./vertical-swap-four-clusters.cdb/data/member_"+str(i)+".vti"
    simplificationThreshold = 0.0
    fieldName = "scalar"
    treeType = int(1)
    tnx = getMergeTree(path,simplificationThreshold,fieldName,treeType)
    
    rootID = -1
    for i in range(len(list(tnx.nodes))):
        if(treeType==0 and tnx.nodes(data=True)[i]["type"]==3):
            rootID = i
        if(treeType==1 and tnx.nodes(data=True)[i]["type"]==0):
            rootID = i
            
    computeBranchDecomposition(tnx, rootID)

    mergeTrees.append((tnx,rootID))

# initialize distance matrices

distanceMatrix_branch = np.zeros((20,20))
distanceMatrix_path = np.zeros((20,20))
distanceMatrix_constrained = np.zeros((20,20))

# compute distances between all pairs of merge trees

for k in range(0,20):
    t1_nx,rootID1 = mergeTrees[k]
    t1_nx_rooted = nx.dfs_tree(t1_nx, rootID1)
    t1_nx_rooted.add_nodes_from((i, t1_nx.nodes[i]) for i in t1_nx_rooted.nodes)
    for l in range(k+1,20):
        t2_nx,rootID2 = mergeTrees[l]
        t2_nx_rooted = nx.dfs_tree(t2_nx, rootID2)
        t2_nx_rooted.add_nodes_from((i, t2_nx.nodes[i]) for i in t2_nx_rooted.nodes)

        print("")
        print("----")
        print(k,l)

        pmd = pathMappingDistance(t1_nx_rooted,rootID1,t2_nx_rooted,rootID2,traceback=False)
        bmd = branchMappingDistance(t1_nx_rooted,rootID1,t2_nx_rooted,rootID2,cost_wasserstein_branch_squared,sqrt=True,traceback=False)
        mted = editDistance_constrained(t1_nx_rooted,rootID1,t2_nx_rooted,rootID2,cost_wasserstein_branch_squared,sqrt=True,traceback=False)
        
        print("Path Mapping Distance:",pmd)
        print("Branch Mapping Distance:",bmd)
        print("Constrained Edit Distance:",mted)
        print("----")
        print("")

        distanceMatrix_path[k,l] = pmd
        distanceMatrix_path[l,k] = pmd
        distanceMatrix_branch[k,l] = bmd
        distanceMatrix_branch[l,k] = bmd
        distanceMatrix_constrained[k,l] = mted
        distanceMatrix_constrained[l,k] = mted

# commpute clustermap for each distance matrix

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_path), method='average')
cm = sns.clustermap(distanceMatrix_path, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_path = cm.dendrogram_row.reordered_ind

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_branch), method='average')
cm = sns.clustermap(distanceMatrix_branch, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_branch = cm.dendrogram_row.reordered_ind

linkage = hc.linkage(sp.distance.squareform(distanceMatrix_constrained), method='average')
cm = sns.clustermap(distanceMatrix_constrained, row_linkage=linkage, col_linkage=linkage, cmap='inferno_r')
reorder_constrained = cm.dendrogram_row.reordered_ind

# save clustermaps as SVG files

fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_path[:, reorder_path][reorder_path], cmap='inferno_r', interpolation='nearest')
plt.axis('off')
plt.colorbar(pos,ax=[ax],location='left',shrink=0.85)
fig.savefig("./clustermap_path_cluster.svg",transparent=True,bbox_inches='tight',pad_inches=0)

fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_branch[:, reorder_branch][reorder_branch], cmap='inferno_r', interpolation='nearest')
plt.axis('off')
plt.colorbar(pos,ax=[ax],location='left',shrink=0.85)
fig.savefig("./clustermap_branch_cluster.svg",transparent=True,bbox_inches='tight',pad_inches=0)

fig, ax = plt.subplots()
pos = ax.imshow(distanceMatrix_constrained[:, reorder_constrained][reorder_constrained], cmap='inferno_r', interpolation='nearest')
plt.axis('off')
plt.colorbar(pos,ax=[ax],location='left',shrink=0.85)
fig.savefig("./clustermap_constrained_cluster.svg",transparent=True,bbox_inches='tight',pad_inches=0)

exit()
