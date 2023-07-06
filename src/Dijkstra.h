#ifndef OCTOMAPDIJKSTRA_DIJKSTRA_H
#define OCTOMAPDIJKSTRA_DIJKSTRA_H

#include <vector>

double Dijkstra(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
         int SID, int FID, std::vector<double> &pathIndexes);

#endif //OCTOMAPDIJKSTRA_DIJKSTRA_H