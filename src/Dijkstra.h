#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <vector>

extern double
Dijkstra(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
         int SID, int FID, std::vector<double> &pathIndexes);

#endif