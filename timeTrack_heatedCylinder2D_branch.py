import numpy as np
import math
    
import matplotlib.pyplot as plt
import vtk

from contourMergeTrees_helpers import *
from baseMetrics import *
from mergeTreeEdit_branch import *

mergeTrees = []
for sim in range(5,6):
    for i in range(0,14):
        t = round(3.0+i*0.10,2)
        ts = "{:.2f}".format(t)
        ts2 = str(int(t*100))
        s = str(int(sim))
        if(sim<10):
            s = "0"+s
        print(s)
        print(ts)

        path = "./heatedCylinder2d.cdb/"+s+"/"+ts+".vti"
        simplificationThreshold = 0.05
        fieldName = "nrrd"
        treeType = int(1)

        nodeScalars,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences,mt_vtk,field = getMergeTree_withVTK(path,simplificationThreshold,fieldName,treeType)
        mergeTrees.append((nodeScalars,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences,mt_vtk,field))

distances = []
matchings = []
for t1 in range(len(mergeTrees)-1):
    t2 = t1+1
    
    nodeScalars1,nodeTypes1,nodePositions1,upIDs1,downIDs1,regionSizes1,persistences1,mt1_vtk,field1 = mergeTrees[t1]
    nodeScalars2,nodeTypes2,nodePositions2,upIDs2,downIDs2,regionSizes2,persistences2,mt2_vtk,field2 = mergeTrees[t2]
        
    rootID = -1
    children1 = []
    for i in range(len(nodeScalars1)):
        if(treeType==0 and nodeTypes1[i]==3):
            rootID1 = i
        if(treeType==1 and nodeTypes1[i]==0):
            rootID1 = i
        children1.append([])
    for i in range(len(upIDs1)):
        downId = downIDs1[i]
        upId = upIDs1[i]
        children1[upId].append(downId)
    for i in range(len(nodeScalars1)):
        print(str(i)+": "+str(children1[i]))
    #print(rootID1)

    rootID2 = -1
    children2 = []
    for i in range(len(nodeScalars2)):
        if(treeType==0 and nodeTypes2[i]==3):
            rootID2 = i
        if(treeType==1 and nodeTypes2[i]==0):
            rootID2 = i
        children2.append([])
    for i in range(len(upIDs2)):
        downId = downIDs2[i]
        upId = upIDs2[i]
        children2[upId].append(downId)
    for i in range(len(nodeScalars2)):
        print(str(i)+": "+str(children2[i]))
    #print(rootID2)

    dist,match = branchMappingDistance(nodeScalars1,children1,rootID1,nodeScalars2,children2,rootID2,cost_wasserstein_branch,True)

    print(t1,t2,dist)
    print(t1,t2,match)
    
    distances.append(dist)
    matchings.append(match)
    
print(distances)
print(matchings)

i=0
j=0
matchLabels = []
for match in matchings:
    matchLabels1 = [-1] * len(mergeTrees[j][0])
    matchLabels2 = [-1] * len(mergeTrees[j+1][0])
    for m in match:
        lastMatchLabels = matchLabels[j-1][1] if j>0 else [-1] * len(nodeScalars1)
        if((not m[0][0]==-1) and (not m[1][0]==-1)):
            #print(m)
            lastMatchLabel = lastMatchLabels[m[0][0]]
            if(lastMatchLabels[m[0][0]]>=0):
                matchLabels1[m[0][0]] = lastMatchLabel
                matchLabels2[m[1][0]] = lastMatchLabel
            else:
                matchLabels1[m[0][0]] = i
                matchLabels2[m[1][0]] = i 
                i+=1
    matchLabels.append((matchLabels1,matchLabels2))
    j+=1

for t in range(len(matchLabels)):
    print(mergeTrees[t][0])
    matchLabels_t = matchLabels[t][0] if t<len(mergeTrees)-1 else matchLabels[t-1][1]
    print(matchLabels_t)
    print("")
#print(matchLabels)

fields = []
multiblock = vtk.vtkMultiBlockDataSet()
multiblock.SetNumberOfBlocks(len(mergeTrees))
for t in range(len(mergeTrees)):
    nodeScalars,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences,mt_vtk,field = mergeTrees[t]

    matchLabels_t = matchLabels[t][0] if t<len(mergeTrees)-1 else matchLabels[t-1][1]

    downIdArray = mt_vtk.GetCellData().GetScalars("downNodeId")
    segIdArray = mt_vtk.GetCellData().GetScalars("SegmentationId")
    segToMatch = dict()
    for i in range(mt_vtk.GetNumberOfCells()):
        downId = downIdArray.GetComponent(i,0)
        segId = segIdArray.GetComponent(i,0)
        label = matchLabels_t[int(downId)]
        segToMatch[segId] = label
    matchingArray = vtk.vtkIntArray()
    matchingArray.SetNumberOfTuples(field.GetNumberOfPoints())
    matchingArray.SetNumberOfComponents(1)
    matchingArray.SetName("MatchingLabel")
    segIdArray = field.GetPointData().GetArray("SegmentationId")
    for i in range(field.GetNumberOfPoints()):
        segId = segIdArray.GetComponent(i,0)
        label = segToMatch[int(segId)]
        matchingArray.SetComponent(i,0,label)
    field.GetPointData().AddArray(matchingArray)
    fields.append(field)
    multiblock.SetBlock(t,field)
    
writer = vtk.vtkXMLMultiBlockDataWriter()
writer.SetFileName("./output/tracking_heatedCylinderMB.vtm")
writer.SetInputData(multiblock)
writer.Write()

