#pragma once

#include <vector>
#include <set>

float editDistance_path(std::vector<float> &nodes1,
                        std::vector<std::vector<int>> &topo1,
                        int rootID1,
                        std::vector<float> &nodes2,
                        std::vector<std::vector<int>> &topo2,
                        int rootID2);
