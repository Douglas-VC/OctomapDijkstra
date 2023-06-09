#include "Dijkstra.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <vector>
#include <algorithm>

namespace coder {
    void find_connections(::coder::array<int, 1U> &nodeIndex, int counter, const coder::array<int, 2U> &E, double b_I) {
        int aux_counter = 0;
        nodeIndex.set_size(counter);
        for (int index = 0; index < counter - 1; index++) {
            if (E[index] == b_I) {
                nodeIndex[aux_counter] = index;
                aux_counter++;
            }
        }
        nodeIndex.set_size(aux_counter);
    }

    void find_not_NaN(::coder::array<int, 1U> &nodeIndex, int numberOfNodes, coder::array<double, 1U> &iTable) {
        int aux_counter = 0;
        nodeIndex.set_size(numberOfNodes);
        for (int index = 0; index < numberOfNodes - 1; index++) {
            if (!rtIsNaN(iTable[index])) {
                nodeIndex[aux_counter] = index;
                aux_counter++;
            }
        }
        nodeIndex.set_size(aux_counter);
    }
}

double Dijkstra(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
                int SID, int FID, std::vector<double> &pathIndexes) {

    coder::array<double, 1U> iTable;
    coder::array<double, 1U> jTable;
    coder::array<double, 1U> minCost;
    coder::array<int, 2U> E;
    coder::array<int, 1U> nodeIndex;
    coder::array<int, 1U> neighbours;
    coder::array<boolean_T, 1U> auxArray2;
    coder::array<boolean_T, 1U> isSettled;
    double b_I;
    int b_J;
    int nx;

    //  Find the minimum costs and paths using Dijkstra's Algorithm
    //  Initializations
    int numberOfNodes = int(AdjacencyList.size());
    iTable.set_size(numberOfNodes);
    jTable.set_size(numberOfNodes);
    minCost.set_size(numberOfNodes);
    isSettled.set_size(numberOfNodes);
    auxArray2.set_size(numberOfNodes);
    for (int index = 0; index < numberOfNodes; index++) {
        iTable[index] = rtNaN;
        minCost[index] = rtInf;
        isSettled[index] = false;
    }

    b_I = SID;
    minCost[SID] = 0.0;
    iTable[SID] = 0.0;
    isSettled[SID] = true;

    std::vector<std::vector<double>> paths(numberOfNodes, std::vector<double>(0));
    paths[SID].push_back(SID);

    //  Execute Dijkstra's Algorithm for this vertex
    while (!isSettled[FID]) {
        //  Update the table
        jTable = iTable;
        iTable[static_cast<int>(b_I)] = rtNaN;

        //  Calculate the costs to the neighbor nodes and record paths
        nx = int(AdjacencyList[static_cast<int>(b_I)].size());
        for (int index = 0; index < nx; index++) {
            b_J = AdjacencyList[static_cast<int>(b_I)][index];
            if (!isSettled[b_J]) {
                float c = CostList[static_cast<int>(b_I)][index];
                if (rtIsNaN(jTable[b_J]) ||
                    (jTable[b_J] > static_cast<float>(jTable[static_cast<int>(b_I)]) + c)) {
                    iTable[b_J] = static_cast<float>(jTable[static_cast<int>(b_I)]) + c;
                    paths[static_cast<int>(b_J)] = (paths[static_cast<int>(b_I)]);
                    paths[static_cast<int>(b_J)].push_back(b_J);
                } else {
                    iTable[b_J] = jTable[b_J];
                }
            }
        }

        //  Find values in the table
        coder::find_not_NaN(nodeIndex, numberOfNodes, iTable);

        //  Settle the minimum value in the table
        nx = nodeIndex.size(0);

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
