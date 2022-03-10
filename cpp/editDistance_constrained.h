#pragma once

#include <vector>
#include <set>

float editDistance_constrained_recursive(   std::vector<float> &nodes1,
                                            std::vector<std::vector<int>> &topo1,
                                            int rootID1,
                                            std::vector<float> &nodes2,
                                            std::vector<std::vector<int>> &topo2,
                                            int rootID2);
float editDistance_constrained( std::vector<float> &nodes1,
                                std::vector<std::vector<int>> &topo1,
                                int rootID1,
                                std::vector<float> &nodes2,
                                std::vector<std::vector<int>> &topo2,
                                int rootID2);
inline float editCost(float v1, float v2);
float editDistance_constrained_tree(int curr1, int curr2,
                                    std::vector<std::vector<float>> &memT,
                                    std::vector<std::vector<float>> &memF,
                                    std::vector<float> &nodes1,
                                    std::vector<std::vector<int>> &topo1,
                                    int rootID1,
                                    std::vector<float> &nodes2,
                                    std::vector<std::vector<int>> &topo2,
                                    int rootID2);
float editDistance_constrained_forest(  int curr1, int curr2,
                                        std::vector<std::vector<float>> &memT,
                                        std::vector<std::vector<float>> &memF,
                                        std::vector<float> &nodes1,
                                        std::vector<std::vector<int>> &topo1,
                                        int rootID1,
                                        std::vector<float> &nodes2,
                                        std::vector<std::vector<int>> &topo2,
                                        int rootID2);