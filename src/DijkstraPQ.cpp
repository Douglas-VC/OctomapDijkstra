#include "DijkstraPQ.h"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <queue>

using std::vector;
using std::pair;

// Dijkstra implementation using Priority Queue
double DijkstraPQ(const vector<vector<int>> &AdjacencyList, const vector<vector<float>> &CostList,
                int SID, int FID, vector<double> &pathIndexes) {

    int size {static_cast<int>(AdjacencyList.size())};

    std::priority_queue<pair<float, int>, vector<pair<float, int>>, std::greater<>> pq;
    vector<float> cost(size, static_cast<float>(INT32_MAX));
    vector<int> parents(size);

    pq.emplace(0.0f, SID);
    cost[SID] = 0.0f;
    parents[SID] = SID;

    double finalCost {};
    while(!pq.empty()) {
        float currCost {pq.top().first};
        int node {pq.top().second};
        pq.pop();

        if(node == FID) {
            finalCost = currCost;
            break;
        }

        for (int i = 0; i < AdjacencyList[node].size(); ++i) {
            int adjNode {AdjacencyList[node][i]};
            float newCost {currCost + CostList[node][i]};
            if(newCost < cost[adjNode]) {
                cost[adjNode] = newCost;
                pq.emplace(newCost, adjNode);
                parents[adjNode] = node;
            }
        }
    }

    int node {FID};
    while(parents[node] != node) {
        pathIndexes.push_back(node);
        node = parents[node];
    }
    pathIndexes.push_back(SID);
    std::reverse(pathIndexes.begin(), pathIndexes.end());

    return finalCost;
}