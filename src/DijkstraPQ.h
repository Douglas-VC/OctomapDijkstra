#ifndef OCTOMAPDIJKSTRA_DIJKSTRAPQ_H
#define OCTOMAPDIJKSTRA_DIJKSTRAPQ_H

#include <vector>

// Dijkstra implementation using a Priority Queue
float DijkstraPQ(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
                int SID, int FID, std::vector<int> &pathIndexes);

#endif //OCTOMAPDIJKSTRA_DIJKSTRAPQ_H
