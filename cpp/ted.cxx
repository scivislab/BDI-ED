// std includes
#include <iostream>
#include <chrono>
#include <stack>

// ttk includes
#include <ttkFTMTree.h>
#include <ttkTopologicalSimplificationByPersistence.h>

// vtk includes
#include <vtkXMLImageDataReader.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>

// file includes
#include "editDistance_constrained.h"
#include "editDistance_branch.h"

using mergeTreeTriple = std::tuple<std::vector<float>,std::vector<std::vector<int>>,int>;

mergeTreeTriple getMergeTree(vtkSmartPointer<vtkUnstructuredGrid> mt_nodes,vtkSmartPointer<vtkUnstructuredGrid> mt, int treeType, bool scalarLabels = false);

int main() {
    int treeType = 1;
    float simplificationThreshold = 0.000002;
    std::string fieldName = "nrrd";

    vtkNew<vtkXMLImageDataReader> reader;
    reader->SetFileName("./../heatedCylinder2d.cdb/01/4.20.vti");

    vtkNew<ttkTopologicalSimplificationByPersistence> simplification{};
    simplification->SetPersistenceThreshold(simplificationThreshold);
    simplification->SetInputArrayToProcess(0,0,0,0,fieldName.c_str());
    simplification->SetInputConnection(reader->GetOutputPort());

    vtkNew<ttkFTMTree> ftm1{};
    ftm1->SetInputConnection(simplification->GetOutputPort());
    ftm1->SetInputArrayToProcess(0,0,0,0,fieldName.c_str());
    ftm1->SetTreeType(treeType);
    ftm1->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mt1 = vtkUnstructuredGrid::SafeDownCast( ftm1->GetOutput(1) );

    reader->SetFileName("./../heatedCylinder2d.cdb/05/4.20.vti");

    vtkNew<ttkFTMTree> ftm2{};
    ftm2->SetInputConnection(simplification->GetOutputPort());
    ftm2->SetInputArrayToProcess(0,0,0,0,fieldName.c_str());
    ftm2->SetTreeType(treeType);
    ftm2->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mt2 = vtkUnstructuredGrid::SafeDownCast( ftm2->GetOutput(1) );

    std::cout << mt1->GetNumberOfPoints() << " " << mt1->GetNumberOfCells() << std::endl;
    std::cout << mt2->GetNumberOfPoints() << " " << mt2->GetNumberOfCells() << std::endl;

    auto t1 = getMergeTree(vtkUnstructuredGrid::SafeDownCast(ftm1->GetOutput(0)),vtkUnstructuredGrid::SafeDownCast(ftm1->GetOutput(1)),treeType);
    auto t2 = getMergeTree(vtkUnstructuredGrid::SafeDownCast(ftm2->GetOutput(0)),vtkUnstructuredGrid::SafeDownCast(ftm2->GetOutput(1)),treeType);

    auto t1_s = getMergeTree(vtkUnstructuredGrid::SafeDownCast(ftm1->GetOutput(0)),vtkUnstructuredGrid::SafeDownCast(ftm1->GetOutput(1)),treeType,true);
    auto t2_s = getMergeTree(vtkUnstructuredGrid::SafeDownCast(ftm2->GetOutput(0)),vtkUnstructuredGrid::SafeDownCast(ftm2->GetOutput(1)),treeType,true);

    auto nodes1 = std::get<0>(t1);
    auto topo1 = std::get<1>(t1);
    auto rootID1 = std::get<2>(t1);
    auto nodes2 = std::get<0>(t2);
    auto topo2 = std::get<1>(t2);
    auto rootID2 = std::get<2>(t2);

    auto nodes1_s = std::get<0>(t1_s);
    auto topo1_s = std::get<1>(t1_s);
    auto rootID1_s = std::get<2>(t1_s);
    auto nodes2_s = std::get<0>(t2_s);
    auto topo2_s = std::get<1>(t2_s);
    auto rootID2_s = std::get<2>(t2_s);

    auto start = std::chrono::high_resolution_clock::now();
    auto dist = editDistance_constrained(nodes1,topo1,rootID1,nodes2,topo2,rootID2);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "d_c: " << dist << std::endl;
    std::cout << "Time d_c: " << duration.count()*0.000001 << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    dist = editDistance_branch(nodes1_s,topo1_s,rootID1_s,nodes2_s,topo2_s,rootID2_s);
    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "d_b bottom-up: " << dist << std::endl;
    std::cout << "Time d_b bottom-up: " << duration.count()*0.000001 << " seconds" << std::endl;

    return 0;
}

mergeTreeTriple getMergeTree(vtkSmartPointer<vtkUnstructuredGrid> mt_nodes,vtkSmartPointer<vtkUnstructuredGrid> mt, int treeType, bool scalarLabels){
    std::map<std::tuple<float,float,float>,float> posToScalarMap;
    std::vector<float> nodeScalars(mt->GetNumberOfPoints());
    std::vector<int> nodeTypes(mt->GetNumberOfPoints());
    std::vector<int> upIDs(mt->GetNumberOfCells());
    std::vector<int> downIDs(mt->GetNumberOfCells());
    std::vector<float> regionSizes(mt->GetNumberOfCells());
    std::vector<float> persistences(mt->GetNumberOfCells());
    std::vector<float> nodePersistences(mt->GetNumberOfPoints(),1000.);
    for (int i=0; i<mt->GetNumberOfPoints(); i++){
        int nid = mt_nodes->GetPointData()->GetArray("NodeId")->GetComponent(i,0);
        int nt = mt_nodes->GetPointData()->GetArray("CriticalType")->GetComponent(i,0);
        float ns = mt_nodes->GetPointData()->GetArray("Scalar")->GetComponent(nid,0);
        nodeTypes[nid] = nt;
        nodeScalars[nid] = ns;
    }
    for (int i=0; i<mt->GetNumberOfCells(); i++){
        int downId = int(mt->GetCellData()->GetScalars("downNodeId")->GetComponent(i,0));
        int upId = int(mt->GetCellData()->GetScalars("upNodeId")->GetComponent(i,0));
        float rs = mt->GetCellData()->GetScalars("RegionSize")->GetComponent(i,0);
        float pers = abs(nodeScalars[upId] - nodeScalars[downId]);
        upIDs[i] = upId;
        downIDs[i] = downId;
        regionSizes[i] = rs;
        persistences[i] = pers;
    }
    
    int rootID = -1;
    std::vector<std::vector<int>> children;
    for (int i=0; i<nodeScalars.size(); i++){
        if(treeType==0 and nodeTypes[i]==3){
            rootID = i;
        }
        if(treeType==1 and nodeTypes[i]==0){
            rootID = i;
        }
        children.push_back(std::vector<int>());
    }
    for (int i=0; i<upIDs.size(); i++){
        int downId = downIDs[i];
        int upId = upIDs[i];
        children[upId].push_back(downId);
        nodePersistences[downId] = persistences[i];//std::abs(nodeScalars[downId]-nodeScalars[upId]);
    }


    // transform to dfs order
    std::vector<int> dfs_reordering(children.size());
    for(int i=0; i<children.size(); i++) dfs_reordering[i] = i;
    std::stack<int> stack;
    stack.push(rootID);
    int currIdx = children.size()-1;
    while(!stack.empty()){
        int nIdx = stack.top();
        stack.pop();
        dfs_reordering[nIdx] = currIdx;
        currIdx--;
        for(int cIdx : children[nIdx]){
            stack.push(cIdx);
        }
    }
    std::vector<std::vector<int>> children_dfs(children.size());
    std::vector<float> labels_dfs(nodePersistences.size());
    std::vector<float> scalars_dfs(nodeScalars.size());
    for(int i=0; i<children.size(); i++){
        int dfsIdx = dfs_reordering[i];
        labels_dfs[dfsIdx] = nodePersistences[i];
        scalars_dfs[dfsIdx] = nodeScalars[i];
        for(int c : children[i]){
            children_dfs[dfsIdx].push_back(dfs_reordering[c]);
        }
    }
    
    return std::make_tuple(scalarLabels ? scalars_dfs : labels_dfs,children_dfs,labels_dfs.size()-1);
}
