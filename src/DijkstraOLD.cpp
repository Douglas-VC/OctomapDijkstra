#include "DijkstraOLD.h"
#include <vector>
#include <algorithm>
#include <cstdint>

using std::vector;

void find_not_NaN(vector<int> &nodeIndex, int numberOfNodes, vector<float> &iTable) {
    int aux_counter = 0;
    nodeIndex.resize(numberOfNodes);
    for (int index = 0; index < numberOfNodes - 1; index++) {
        if (iTable[index] != -1) {
            nodeIndex[aux_counter] = index;
            aux_counter++;
        }
    }
    nodeIndex.resize(aux_counter);
}

// OLD Dijkstra implementation
float DijkstraOLD(const vector<vector<int>> &AdjacencyList, const vector<vector<float>> &CostList,
                int SID, int FID, vector<int> &pathIndexes) {

    int numberOfNodes {static_cast<int>(AdjacencyList.size())};
    vector<float> iTable(numberOfNodes, -1);
    vector<float> jTable(numberOfNodes, -1);
    vector<float> minCost(numberOfNodes, static_cast<float>(INT32_MAX));
    vector<int> nodeIndex(numberOfNodes);
    vector<int> neighbours(numberOfNodes);
    vector<bool> isSettled(numberOfNodes);
    vector<vector<int>> E;

    int b_I;
    int b_J;
    int nx;

    b_I = SID;
    minCost[SID] = 0.0;
    iTable[SID] = 0.0;
    isSettled[SID] = true;

    vector<vector<int>> paths(numberOfNodes, vector<int>());
    paths[SID].push_back(SID);

    //  Execute Dijkstra's Algorithm for this vertex
    while (!isSettled[FID]) {
        //  Update the table
        jTable = iTable;
        iTable[b_I] = -1;

        //  Calculate the costs to the neighbor nodes and record paths
        nx = int(AdjacencyList[b_I].size());
        for (int index = 0; index < nx; index++) {
            b_J = AdjacencyList[b_I][index];
            if (!isSettled[b_J]) {
                float c = CostList[b_I][index];
                if (jTable[b_J] == -1 || jTable[b_J] > jTable[b_I] + c) {
                    iTable[b_J] = jTable[b_I] + c;
                    paths[b_J] = (paths[b_I]);
                    paths[b_J].push_back(b_J);
                } else {
                    iTable[b_J] = jTable[b_J];
                }
            }
        }
        //  Find values in the table
        find_not_NaN(nodeIndex, numberOfNodes, iTable);

        //  Settle the minimum value in the table
        nx = static_cast<int>(nodeIndex.size());

        float minC = iTable[nodeIndex[0]];

        int aux = 0;
        float d;
        for (int index = 1; index < nx; index++) {
            d = iTable[nodeIndex[index]];
            if (minC > d) {
                minC = d;
                aux = index;
            }
        }

        b_I = nodeIndex[aux];
        minCost[b_I] = iTable[b_I];
        isSettled[b_I] = true;
    }

    //  Store costs and paths
    pathIndexes = paths[FID];
    return minCost[FID];
}
