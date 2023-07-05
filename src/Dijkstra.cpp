#include "Dijkstra.h"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <iostream>

using std::vector;

//void find_connections(vector<int> &nodeIndex, int counter, vector<vector<int>> &E, double b_I) {
//    int aux_counter = 0;
//    nodeIndex.resize(counter);
//    for (int index = 0; index < counter - 1; index++) {
//        if (E[index] == static_cast<int>(b_I)) {
//            nodeIndex[aux_counter] = index;
//            aux_counter++;
//        }
//    }
//    nodeIndex.set_size(aux_counter);
//}

void find_not_NaN(vector<int> &nodeIndex, int numberOfNodes, vector<double> &iTable) {
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

double Dijkstra(const vector<vector<int>> &AdjacencyList, const vector<vector<float>> &CostList,
                int SID, int FID, vector<double> &pathIndexes) {

    int numberOfNodes {static_cast<int>(AdjacencyList.size())};
    vector<double> iTable(numberOfNodes, -1);
    vector<double> jTable(numberOfNodes, -1);
    vector<double> minCost(numberOfNodes, INT32_MAX);
    vector<int> nodeIndex(numberOfNodes);
    vector<int> neighbours(numberOfNodes);
    vector<bool> isSettled(numberOfNodes);
    vector<vector<int>> E;

    double b_I;
    int b_J;
    int nx;

    b_I = SID;
    minCost[SID] = 0.0;
    iTable[SID] = 0.0;
    isSettled[SID] = true;

    vector<vector<double>> paths(numberOfNodes, vector<double>());
    paths[SID].push_back(SID);

    //  Execute Dijkstra's Algorithm for this vertex
    while (!isSettled[FID]) {
        //  Update the table
        jTable = iTable;
        iTable[static_cast<int>(b_I)] = -1;

        //  Calculate the costs to the neighbor nodes and record paths
        nx = int(AdjacencyList[static_cast<int>(b_I)].size());
        for (int index = 0; index < nx; index++) {
            b_J = AdjacencyList[static_cast<int>(b_I)][index];
            if (!isSettled[b_J]) {
                float c = CostList[static_cast<int>(b_I)][index];
                if (jTable[b_J] == -1 || jTable[b_J] > static_cast<float>(jTable[static_cast<int>(b_I)]) + c) {
                    iTable[b_J] = static_cast<float>(jTable[static_cast<int>(b_I)]) + c;
                    paths[static_cast<int>(b_J)] = (paths[static_cast<int>(b_I)]);
                    paths[static_cast<int>(b_J)].push_back(b_J);
                } else {
                    iTable[b_J] = jTable[b_J];
                }
            }
        }
        //  Find values in the table
        find_not_NaN(nodeIndex, numberOfNodes, iTable);

        //  Settle the minimum value in the table
        nx = nodeIndex.size();

        b_I = iTable[static_cast<int>(nodeIndex[0])];

        int aux = 0;
        double d;
        for (int index = 1; index < nx; index++) {
            d = iTable[static_cast<int>(nodeIndex[index])];
            if (b_I > d) {
                b_I = d;
                aux = index;
            }
        }

        b_I = nodeIndex[aux];
        minCost[static_cast<int>(b_I)] = iTable[static_cast<int>(b_I)];
        isSettled[static_cast<int>(b_I)] = true;
    }

    //  Store costs and paths
    pathIndexes = paths[FID];
    return minCost[FID];
}
