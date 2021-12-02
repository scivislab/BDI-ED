import numpy as np
import math

import vtk
import topologytoolkit as ttk
    
import matplotlib.pyplot as plt
import networkx as nx


def getMergeTree(path,simplificationThreshold,fieldName,treeType):
    reader = vtk.vtkXMLImageDataReader()
    #reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)

    # reader.Update()
    # test = reader.GetOutput()
    # print(test)

    simplification = ttk.ttkTopologicalSimplificationByPersistence()
    simplification.SetPersistenceThreshold(simplificationThreshold)
    simplification.SetInputArrayToProcess(0,0,0,0,fieldName)
    simplification.SetInputConnection(reader.GetOutputPort())
    simplification.SetDebugLevel(3)

    ftm = ttk.ttkFTMTree()
    ftm.SetDebugLevel(3)
    ftm.SetInputArrayToProcess(0,0,0,0,fieldName)
    ftm.SetInputConnection(simplification.GetOutputPort())
    ftm.SetTreeType(treeType)

    ftm.Update()
    mt = ftm.GetOutput(1)
    mt_nodes = ftm.GetOutput(0)
    posToScalarMap = dict()
    nodeScalars = []
    nodeTypes = []
    upIDs = []
    downIDs = []
    regionSizes = []
    persistences = []
    nodePositions = []
    for i in range(mt.GetNumberOfPoints()):
        pos1 = mt.GetPoint(i)
        posToScalarMap[pos1] = mt.GetPointData().GetArray("Scalar").GetComponent(i,0)
        nodeScalars.append(0)
        nodeTypes.append(0)
        nodePositions.append(pos1)
    for i in range(mt.GetNumberOfPoints()):
        nid = int(mt_nodes.GetPointData().GetArray("NodeId").GetComponent(i,0))
        nt = int(mt_nodes.GetPointData().GetArray("CriticalType").GetComponent(i,0))
        pos = mt_nodes.GetPoint(i)
        ns = posToScalarMap[pos]
        nodeTypes[nid] = nt
        nodeScalars[nid] = ns
        nodePositions[nid] = pos
    for i in range(mt.GetNumberOfCells()):
        downId = int(mt.GetCellData().GetScalars("downNodeId").GetComponent(i,0))
        upId = int(mt.GetCellData().GetScalars("upNodeId").GetComponent(i,0))
        rs = mt.GetCellData().GetScalars("RegionSize").GetComponent(i,0)
        pers = abs(nodeScalars[upId] - nodeScalars[downId])
        upIDs.append(upId)
        downIDs.append(downId)
        regionSizes.append(rs)
        persistences.append(pers)
    
    return nodeScalars,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences
        
def getMergeTree_withVTK(path,simplificationThreshold,fieldName,treeType):
    reader = vtk.vtkXMLImageDataReader()
    #reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)

    # reader.Update()
    # test = reader.GetOutput()
    # print(test)

    simplification = ttk.ttkTopologicalSimplificationByPersistence()
    simplification.SetPersistenceThreshold(simplificationThreshold)
    simplification.SetInputArrayToProcess(0,0,0,0,fieldName)
    simplification.SetInputConnection(reader.GetOutputPort())
    simplification.SetDebugLevel(3)

    ftm = ttk.ttkFTMTree()
    ftm.SetDebugLevel(3)
    ftm.SetInputArrayToProcess(0,0,0,0,fieldName)
    ftm.SetInputConnection(simplification.GetOutputPort())
    ftm.SetTreeType(treeType)

    ftm.Update()
    mt = ftm.GetOutput(1)
    mt_nodes = ftm.GetOutput(0)
    posToScalarMap = dict()
    nodeScalars = []
    nodeTypes = []
    upIDs = []
    downIDs = []
    regionSizes = []
    persistences = []
    nodePositions = []
    for i in range(mt.GetNumberOfPoints()):
        pos1 = mt.GetPoint(i)
        posToScalarMap[pos1] = mt.GetPointData().GetArray("Scalar").GetComponent(i,0)
        nodeScalars.append(0)
        nodeTypes.append(0)
        nodePositions.append(pos1)
    for i in range(mt.GetNumberOfPoints()):
        nid = int(mt_nodes.GetPointData().GetArray("NodeId").GetComponent(i,0))
        nt = int(mt_nodes.GetPointData().GetArray("CriticalType").GetComponent(i,0))
        pos = mt_nodes.GetPoint(i)
        ns = posToScalarMap[pos]
        nodeTypes[nid] = nt
        nodeScalars[nid] = ns
        nodePositions[nid] = pos
    for i in range(mt.GetNumberOfCells()):
        downId = int(mt.GetCellData().GetScalars("downNodeId").GetComponent(i,0))
        upId = int(mt.GetCellData().GetScalars("upNodeId").GetComponent(i,0))
        rs = mt.GetCellData().GetScalars("RegionSize").GetComponent(i,0)
        pers = abs(nodeScalars[upId] - nodeScalars[downId])
        upIDs.append(upId)
        downIDs.append(downId)
        regionSizes.append(rs)
        persistences.append(pers)
    
    return nodeScalars,nodeTypes,nodePositions,upIDs,downIDs,regionSizes,persistences,mt,ftm.GetOutput(2)

def getAllMatchings(l1,l2):
    if(len(l1)==0 and len(l2)==0):
        return set([(frozenset([]),frozenset([]),frozenset([]))])
    else:
        m = set()
        for e1 in l1:
            for e2 in l2:
                l1_ = l1.copy()
                l2_ = l2.copy()
                l1_.remove(e1)
                l2_.remove(e2)
                m_ = getAllMatchings(l1_,l2_)
                m__ = []
                for x in m_:
                    x_ = list(x[1])
                    x_.append((e1,e2))
                    x__ = frozenset(x_)
                    m__.append((x[0],x__,x[2]))
                m = m.union(set(m__))
        for e1 in l1:
            l1_ = l1.copy()
            l1_.remove(e1)
            m_ = getAllMatchings(l1_,l2)
            m__ = []
            for x in m_:
                x_ = list(x[0])
                x_.append((e1))
                x__ = frozenset(x_)
                m__.append((x__,x[1],x[2]))
            m = m.union(set(m__))
        for e2 in l2:
            l2_ = l2.copy()
            l2_.remove(e2)
            m_ = getAllMatchings(l1,l2_)
            m__ = []
            for x in m_:
                x_ = list(x[2])
                x_.append((e2))
                x__ = frozenset(x_)
                m__.append((x[0],x[1],x__))
            m = m.union(set(m__))
        return m

def getBranchDecomposition(topo,nodeScalars,rootIdx):
    nodeBranches = list(range(len(nodeScalars)))
    branchLabels = dict()
    computeBranches(topo,nodeScalars,nodeBranches,branchLabels,rootIdx)
    return nodeBranches,branchLabels

def computeBranches(topo,nodeScalars,nodeBranches,branchLabels,curr):
    if(len(topo[curr])==0):
        nodeBranches[curr] = curr
        return
    optPers = -1.0
    for child in topo[curr]:
        computeBranches(topo,nodeScalars,nodeBranches,branchLabels,child)
        pers = abs(nodeScalars[nodeBranches[child]]-nodeScalars[curr])
        if(pers > optPers):
            optPers = pers
            nodeBranches[curr] = nodeBranches[child]
        branchLabels[nodeBranches[child]] = (nodeScalars[nodeBranches[child]],nodeScalars[curr])
