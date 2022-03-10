#include "editDistance_constrained.h"
#include "munkres.hpp"

#include <iostream>
#include <vector>
#include <tuple>
#include <limits>
#include <cmath>
#include <float.h>
#include <set>
#include <algorithm>

float editDistance_constrained( std::vector<float> &nodes1,
                                std::vector<std::vector<int>> &topo1,
                                int rootID1,
                                std::vector<float> &nodes2,
                                std::vector<std::vector<int>> &topo2,
                                int rootID2){

    // initialize memoization tables
    int nn1 = nodes1.size();
    int nn2 = nodes2.size();
    int dim1 = nn1+1;
    int dim2 = nn2+1;
    
    std::vector<float> memT(dim1*dim2,-1);
    std::vector<float> memF(dim1*dim2,-1);
    //float* memT = new float[dim1*dim2];
    //float* memF = new float[dim1*dim2];

    float nan = std::numeric_limits<float>::quiet_NaN();

    std::vector<float> memEdit(dim1*dim2,-1);
    //float* memEdit = new float[dim1*dim2];
    memEdit[nn1+nn2*dim1] = 0;
    for(int i=0; i<nn1; i++){
        memEdit[i+nn2*dim1] = editCost(nodes1[i],nan);
    }
    for(int j=0; j<nn2; j++){
        memEdit[nn1+j*dim1] = editCost(nan,nodes2[j]);
    }
    for(int i=0; i<nn1; i++){
        for(int j=0; j<nn2; j++){
            memEdit[i+j*dim1] = editCost(nodes1[i],nodes2[j]);
        }
    }

    memT[nn1+nn2*dim1] = 0;
    memF[nn1+nn2*dim1] = 0;
    for(int i=0; i<nn1; i++){
        memF[i+nn2*dim1] = 0;
        for(int c : topo1[i]){
            memF[i+nn2*dim1] += memT[c+nn2*dim1];
        }
        memT[i+nn2*dim1] = memF[i+nn2*dim1] + memEdit[i+nn2*dim1];
    }
    for(int j=0; j<nn2; j++){
        memF[nn1+j*dim1] = 0;
        for(int c : topo2[j]){
            memF[nn1+j*dim1] += memT[nn1+c*dim1];
        }
        memT[nn1+j*dim1] = memF[nn1+j*dim1] + memEdit[nn1+j*dim1];
    }
    for(int curr1=0; curr1<nn1; curr1++){
        for(int curr2=0; curr2<nn2; curr2++){
            {
                float d = FLT_MAX;
                int cn1 = topo1[curr1].size();
                int cn2 = topo2[curr2].size();
                // Cases for mapping some trees in first forest to some trees in second forest and deleting all other trees in both forests
                if(cn1<=2 && cn2<=2){
                    int c11 = cn1>0 ? topo1[curr1][0] : nn1;
                    int c12 = cn1>1 ? topo1[curr1][1] : nn1;
                    int c21 = cn2>0 ? topo2[curr2][0] : nn2;
                    int c22 = cn2>1 ? topo2[curr2][1] : nn2;
                    d = std::min(d,memT[c11+c21*dim1]+memT[c12+c22*dim1]);
                    d = std::min(d,memT[c11+c22*dim1]+memT[c12+c21*dim1]);
                }
                else{
                    auto f = [&] (unsigned r, unsigned c) {
                        int c1 = r<topo1[curr1].size() ? topo1[curr1][r] : nn1;
                        int c2 = c<topo2[curr2].size() ? topo2[curr2][c] : nn2;
                        return memT[c1+c2*dim1];
                    };
                    int size = std::max(topo1[curr1].size(),topo2[curr2].size()) + 1;
                    auto matching = munkres_algorithm<float>(size, size, f);
                    float d_ = 0;
                    for(auto m : matching) d_ += f(m.first,m.second);
                    d = std::min(d,d_);
                }
                // Cases for mapping one subtree in first forest to second forest and deleting all other trees in first forest 
                for (int child1 : topo1[curr1]){
                    if(topo1[child1].size()==0){
                        continue;
                    }
                    float d_ = (memF[curr1+nn2*dim1]
                                + memF[child1+curr2*dim1]
                                - memF[child1+nn2*dim1]);
                    d = std::min(d,d_);
                }
                // Cases for mapping one subtree in second forest to first forest and deleting all other trees in second forest 
                for (int child2 : topo2[curr2]){
                    if(topo2[child2].size()==0){
                        continue;
                    }
                    float d_ = (memF[nn1+curr2*dim1]
                                + memF[curr1+child2*dim1]
                                - memF[nn1+child2*dim1]);
                    d = std::min(d,d_);
                }
                memF[curr1+curr2*dim1] = d;
            }
            {
                // Case for matching curr1 to curr2 and recurse to subforests
                float d = (memF[curr1+curr2*dim1]
                            + memEdit[curr1+curr2*dim1]);
                // Cases for deleting curr1 and mapping one child to curr2
                for (int child1 : topo1[curr1]){
                    float d_ = (memT[curr1+nn2*dim1]
                                + memT[child1+curr2*dim1]
                                - memT[child1+nn2*dim1]);
                    d = std::min(d,d_);
                }
                // Cases for deleting curr2 and mapping one child to curr1
                for (int child2 : topo2[curr2]){
                    float d_ = (memT[nn1+curr2*dim1]
                                + memT[curr1+child2*dim1]
                                - memT[nn1+child2*dim1]);
                    d = std::min(d,d_);
                }
                memT[curr1+curr2*dim1] = d;
            }
        }
    }

    float res = memT[rootID1+rootID2*dim1];
    // delete[] memT;
    // delete[] memF;
    // delete[] memEdit;
    return res;
}

float editDistance_constrained_recursive(   std::vector<float> &nodes1,
                                            std::vector<std::vector<int>> &topo1,
                                            int rootID1,
                                            std::vector<float> &nodes2,
                                            std::vector<std::vector<int>> &topo2,
                                            int rootID2){

    // initialize memoization tables
    std::vector<std::vector<float>> memT(nodes1.size()+1,std::vector<float>(nodes2.size()+1,-1));
    std::vector<std::vector<float>> memF(nodes1.size()+1,std::vector<float>(nodes2.size()+1,-1));
    //std::map<std::pair<int,int>,float> memT;
    //std::map<std::pair<int,int>,float> memF;

    float d = editDistance_constrained_tree(rootID1,rootID2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2);

    return d;
}

inline float editCost(float v1, float v2){
    if(std::isnan(v1)) return v2;
    if(std::isnan(v2)) return v1;
    return std::abs(v1-v2);
}

// ===================================================================
// Recursive helper function that computes edit distance between two subtrees rooted in curr1,curr2
float editDistance_constrained_tree(int curr1, int curr2,
                                    std::vector<std::vector<float>> &memT,
                                    std::vector<std::vector<float>> &memF,
                                    std::vector<float> &nodes1,
                                    std::vector<std::vector<int>> &topo1,
                                    int rootID1,
                                    std::vector<float> &nodes2,
                                    std::vector<std::vector<int>> &topo2,
                                    int rootID2){

    // If both trees empty, return distance 0
    if(curr1<0 and curr2<0){
        return 0;
    }
    // If first tree empty, return deletion cost of curr2 plus deletion cost of its subforest
    if(curr1<0){
        if(memT[nodes1.size()][curr2]<0){
            float d = (editDistance_constrained_forest(curr1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                     + editCost(std::nanf(""),nodes2[curr2]));
            memT[nodes1.size()][curr2] = d;
        }
        return memT[nodes1.size()][curr2];
    }
    // If second tree empty, return deletion cost of curr1 plus deletion cost of its subforest
    if(curr2<0){
        if(memT[curr1][nodes2.size()]<0){
            float d = (editDistance_constrained_forest(curr1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                     + editCost(nodes1[curr1],std::nanf("")));
            memT[curr1][nodes2.size()] = d;
        }
        return memT[curr1][nodes2.size()];
    }
    // If both trees not empty, find optimal edit operation
    if(memT[curr1][curr2]<0){
        // Case for matching curr1 to curr2 and recurse to subforests
        float d = (editDistance_constrained_forest(curr1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                 + editCost(nodes1[curr1],nodes2[curr2]));
        // Cases for deleting curr1 and mapping one child to curr2
        for (int child1 : topo1[curr1]){
            float d_ = (editDistance_constrained_tree(curr1,-1,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      + editDistance_constrained_tree(child1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      - editDistance_constrained_tree(child1,-1,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2));
            d = std::min(d,d_);
        }
        // Cases for deleting curr2 and mapping one child to curr1
        for (int child2 : topo2[curr2]){
            float d_ = (editDistance_constrained_tree(-1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      + editDistance_constrained_tree(curr1,child2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      - editDistance_constrained_tree(-1,child2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2));
            d = std::min(d,d_);
        }
        memT[curr1][curr2] = d;
    }
    return memT[curr1][curr2];
};

// ===================================================================
// Recursive helper function that computes edit distance between two subforests rooted in curr1,curr2
float editDistance_constrained_forest(  int curr1, int curr2,
                                        std::vector<std::vector<float>> &memT,
                                        std::vector<std::vector<float>> &memF,
                                        std::vector<float> &nodes1,
                                        std::vector<std::vector<int>> &topo1,
                                        int rootID1,
                                        std::vector<float> &nodes2,
                                        std::vector<std::vector<int>> &topo2,
                                        int rootID2){

    // If both trees empty, return distance 0
    if(curr1<0 and curr2<0){
        return 0;
    }
    // If first forest empty, return deletion cost of all trees in second forest
    if(curr1<0){
        if(memF[nodes1.size()][curr2]<0){
            float d = 0.0;
            for (int child2 : topo2[curr2]){
                d += editDistance_constrained_tree(curr1,child2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2);
            }
            memF[nodes1.size()][curr2] = d;
        }
        return memF[nodes1.size()][curr2];
    }
    // If second forest empty, return deletion cost of all trees in first forest
    if(curr2<0){
        if(memF[curr1][nodes2.size()]<0){
            float d = 0.0;
            for (int child1 : topo1[curr1]){
                d += editDistance_constrained_tree(child1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2);
            }
            memF[curr1][nodes2.size()] = d;
        }
        return memF[curr1][nodes2.size()];
    }
    // If both forests not empty, find optimal edit operations
    if(memF[curr1][curr2]<0){
        float d = FLT_MAX;
        // Cases for mapping some trees in first forest to some trees in second forest and deleting all other trees in both forests
        if(topo1[curr1].size()<=2 && topo2[curr2].size()<=2){
            int c11 = topo1[curr1].size()>0 ? topo1[curr1][0] : -1;
            int c12 = topo1[curr1].size()>1 ? topo1[curr1][1] : -1;
            int c21 = topo2[curr2].size()>0 ? topo2[curr2][0] : -1;
            int c22 = topo2[curr2].size()>1 ? topo2[curr2][1] : -1;
            d = std::min(d,editDistance_constrained_tree(c11,c21,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                          +editDistance_constrained_tree(c12,c22,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2));
            d = std::min(d,editDistance_constrained_tree(c11,c22,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                          +editDistance_constrained_tree(c12,c21,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2));
        }
        else{
            auto f = [&] (unsigned r, unsigned c) {
                int c1 = r<topo1[curr1].size() ? topo1[curr1][r] : -1;
                int c2 = c<topo2[curr2].size() ? topo2[curr2][c] : -1;
                return editDistance_constrained_tree(c1,c2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2);
            };
            int size = std::max(topo1[curr1].size(),topo2[curr2].size()) + 1;
            auto matching = munkres_algorithm<float>(size, size, f);
            float d_ = 0;
            for(auto m : matching) d_ += f(m.first,m.second);
            d = std::min(d,d_);
        }
        // Cases for mapping one subtree in first forest to second forest and deleting all other trees in first forest 
        for (int child1 : topo1[curr1]){
            if(topo1[child1].size()==0){
                continue;
            }
            float d_ = (editDistance_constrained_forest(curr1,-1,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      + editDistance_constrained_forest(child1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      - editDistance_constrained_forest(child1,-1,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2));
            d = std::min(d,d_);
        }
        // Cases for mapping one subtree in second forest to first forest and deleting all other trees in second forest 
        for (int child2 : topo2[curr2]){
            if(topo2[child2].size()==0){
                continue;
            }
            float d_ = (editDistance_constrained_forest(-1,curr2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      + editDistance_constrained_forest(curr1,child2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2)
                      - editDistance_constrained_forest(-1,child2,memT,memF,nodes1,topo1,rootID1,nodes2,topo2,rootID2));
            d = std::min(d,d_);
        }
        memF[curr1][curr2] = d;
    }
    return memF[curr1][curr2];
}