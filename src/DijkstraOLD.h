#ifndef OCTOMAPDIJKSTRA_DIJKSTRAOLD_H
#define OCTOMAPDIJKSTRA_DIJKSTRAOLD_H

#include <vector>

// OLD Dijkstra implementation
float DijkstraOLD(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
         int SID, int FID, std::vector<int> &pathIndexes);

#endif //OCTOMAPDIJKSTRA_DIJKSTRAOLD_H