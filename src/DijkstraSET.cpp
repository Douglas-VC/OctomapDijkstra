#include "DijkstraSET.h"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <set>

using std::vector;
using std::pair;

// Dijkstra implementation using a Set
float DijkstraSET(const vector<vector<int>> &AdjacencyList, const vector<vector<float>> &CostList,
                 int start, int end, vector<int> &pathIndexes) {

    int size {static_cast<int>(AdjacencyList.size())};

    std::set<pair<float, int>> set;
    vector<float> cost(size, static_cast<float>(INT32_MAX));
    vector<int> parents(size);

    set.emplace(0.0f, start);
    cost[start] = 0.0f;
    parents[start] = start;

    float finalCost {};
    while(!set.empty()) {
        float currCost {set.begin()->first};
        int node {set.begin()->second};
        set.erase(*set.begin());

        if(node == end) {
            finalCost = currCost;
            break;
        }

        for (int i = 0; i < AdjacencyList[node].size(); ++i) {
            int adjNode {AdjacencyList[node][i]};
            float newCost {currCost + CostList[node][i]};
            if(newCost < cost[adjNode]) {
                if(cost[adjNode] != static_cast<float>(INT32_MAX))
                    set.erase({cost[adjNode], adjNode});
                cost[adjNode] = newCost;
                set.emplace(newCost, adjNode);
                parents[adjNode] = node;
            }
        }
    }

    int node {end};
    while(parents[node] != node) {
        pathIndexes.push_back(node);
        node = parents[node];
    }
    pathIndexes.push_back(start);
    std::reverse(pathIndexes.begin(), pathIndexes.end());

    return finalCost;
}
