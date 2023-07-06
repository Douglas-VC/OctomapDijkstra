#ifndef OCTOMAPDIJKSTRA_DIJKSTRAPQ_H
#define OCTOMAPDIJKSTRA_DIJKSTRAPQ_H

#include <vector>

// Dijkstra implementation using Priority Queue
double DijkstraPQ(const std::vector<std::vector<int>> &AdjacencyList, const std::vector<std::vector<float>> &CostList,
                int SID, int FID, std::vector<double> &pathIndexes);

#endif //OCTOMAPDIJKSTRA_DIJKSTRAPQ_H
