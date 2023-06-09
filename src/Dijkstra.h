#ifndef DIJKSTRA_H
#define DIJKSTRA_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern double
Dijkstra(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
         int SID, int FID, std::vector<double> &pathIndexes);

#endif