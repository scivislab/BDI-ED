import networkx as nx
import vtk
import topologytoolkit as ttk

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

    ftm = ttk.ttkMergeTree()
    ftm.SetDebugLevel(3)
    ftm.SetInputArrayToProcess(0,0,0,0,fieldName)
    ftm.SetInputConnection(simplification.GetOutputPort())
    ftm.SetTreeType(treeType)

    ftm.Update()
    mt = ftm.GetOutput(1)
    mt_nodes = ftm.GetOutput(0)
    mtnx = nx.Graph()
    for i in range(mt.GetNumberOfPoints()):
        # nid = int(mt_nodes.GetPointData().GetArray("NodeId").GetComponent(i,0))
        # print(nid)
        nt = int(mt_nodes.GetPointData().GetArray("CriticalType").GetComponent(i,0))
        pos = mt_nodes.GetPoint(i)
        ns = mt_nodes.GetPointData().GetArray("Scalar").GetComponent(i,0)
        mtnx.add_node(i,scalar=ns,position=pos,type=nt,birth=0,death=0)
    for i in range(mt.GetNumberOfCells()):
        downId = int(mt.GetCellData().GetScalars("downNodeId").GetComponent(i,0))
        upId = int(mt.GetCellData().GetScalars("upNodeId").GetComponent(i,0))
        rs = mt.GetCellData().GetScalars("RegionSize").GetComponent(i,0)
        pers = abs(mtnx.nodes[upId]['scalar'] - mtnx.nodes[downId]['scalar'])
        mtnx.add_edge(upId,downId,persistence=pers,regionSize=rs)
    
    return mtnx

def computeBranchDecomposition(tree,rootIdx):
    rt = nx.dfs_tree(tree, rootIdx)
    computeBranches(tree,rt,rootIdx)
    return

def computeBranches(t,rt,curr):
    if(len(list(rt.neighbors(curr)))==0):
        return curr
    optPers = -1.0
    furthest = None
    for child in rt.neighbors(curr):
        leaf = computeBranches(t,rt,child)
        pers = abs(t.nodes(data=True)[leaf]['scalar']-t.nodes(data=True)[curr]['scalar'])
        if(pers > optPers):
            optPers = pers
            furthest = leaf
        nx.set_node_attributes(t, {curr:t.nodes(data=True)[leaf]['scalar']}, 'birth')
        nx.set_node_attributes(t, {curr:t.nodes(data=True)[curr]['scalar']}, 'death')
        nx.set_node_attributes(t, {leaf:t.nodes(data=True)[leaf]['scalar']}, 'birth')
        nx.set_node_attributes(t, {leaf:t.nodes(data=True)[curr]['scalar']}, 'death')
    return furthest
