#ifndef OCTOMAPDIJKSTRA_DIJKSTRASET_H
#define OCTOMAPDIJKSTRA_DIJKSTRASET_H

#include <vector>

// Dijkstra implementation using a Set
float DijkstraSET(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
                 int SID, int FID, std::vector<int> &pathIndexes);

#endif //OCTOMAPDIJKSTRA_DIJKSTRASET_H
