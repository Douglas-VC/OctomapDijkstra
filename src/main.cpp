#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/common/common.h>
#include <pcl/octree/octree_search.h>
#include <cmath>
#include <octomap/octomap.h>
#include <chrono>
#include "Dijkstra.h"
#include "DijkstraPQ.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using namespace std::chrono;

using PointCloud = pcl::PointCloud<pcl::PointXYZ>::Ptr;
using Octree = pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>;
using timePoint = std::chrono::high_resolution_clock::time_point;

typedef std::vector<std::vector<int>> AdjacencyList;
typedef std::vector<std::vector<float>> CostList;

void printExecutionTime(timePoint start, timePoint stop, const std::string &mode) {
    duration<double, std::milli> duration = stop - start;
    if (mode == "setup") {
        std::cout << "Setup levou "
                  << duration.count() / 1000
                  << " segundos para executar."
                  << std::endl;
    } else if (mode == "dijkstra") {
        std::cout << "Dijkstra levou "
                  << duration.count() / 1000
                  << " segundos para executar."
                  << std::endl;
    }
}

PointCloud importPointCloud() {
    PointCloud cloud(new pcl::PointCloud<pcl::PointXYZ>);

    if (pcl::io::loadPCDFile<pcl::PointXYZ>(
            "/home/douglas/catkin_ws/src/dijkstra_planning/resources/PointClouds/campinhoremaster3.pcd", *cloud) ==
        -1) {
        PCL_ERROR("Couldn't read point cloud file\n");
        exit(-1);
    }

    std::cout << "Nuvem de pontos carregada: "
              << cloud->width * cloud->height
              << " pontos foram obtidos. "
              << std::endl;

    return cloud;
}

PointCloud importOctomapCreatePointCloud(const std::string &octomapFilePath, float &graphAndOctreeResolution) {
    std::cout << "Importando e processando octomap." << std::endl;

    octomap::AbstractOcTree *otTree = octomap::AbstractOcTree::read(octomapFilePath);
    auto octomap = dynamic_cast<octomap::OcTree *>(otTree);

    graphAndOctreeResolution = static_cast<float>(octomap->getResolution());

    PointCloud cloud(new pcl::PointCloud<pcl::PointXYZ>);

    std::ofstream MyFile("/home/douglas/Área de Trabalho/Dijkstra/frontierOctomap.txt");

    for (octomap::OcTree::leaf_iterator it = octomap->begin_leafs(), end = octomap->end_leafs(); it != end; ++it) {
        if (octomap->isNodeOccupied(*it)) {
            octomap::OcTreeKey neighbour = it.getKey();
            neighbour[2] += 1;
            octomap::OcTreeNode *result = octomap->search(neighbour);
            if (result != nullptr && !octomap->isNodeOccupied(*result)) {
                cloud->push_back(pcl::PointXYZ(float(it.getX()), float(it.getY()), float(it.getZ())));
                MyFile << it.getX() << " "
                       << it.getY() << " "
                       << it.getZ() << std::endl;
            }
        }
    }

    delete otTree;
    return cloud;
}

Octree transformPointCloudToOctree(const PointCloud &cloud, float resolution) {
    Octree octree(resolution);

    octree.setInputCloud(cloud);
    octree.addPointsFromInputCloud();
    return octree;
}

void getMatrixDimensions(const PointCloud &cloud, int &rows, int &columns, float graphResolution) {
    pcl::PointXYZ minPt, maxPt;
    pcl::getMinMax3D(*cloud, minPt, maxPt);

    rows = static_cast<int>(floor(fabs((maxPt.y - minPt.y)) / graphResolution)) + 1;
    columns = static_cast<int>(floor(fabs((maxPt.x - minPt.x)) / graphResolution)) + 1;
}

void populateAdjacencyList(AdjacencyList &AdjacencyList, int rows, int columns) {
    int N = rows * columns;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (j > i + 2 * columns) {
                break;
            } else if ((j == i + 1 && (i + 1) % columns != 0)
                       || ((j == i + columns) && i < N - columns)
                       || (j == i + columns - 1 && i < N - columns && (j + 1) % columns != 0)
                       || (j == i + columns + 1 && i < N - columns && (i + 1) % columns != 0)) {
                AdjacencyList[i].push_back(j);
                AdjacencyList[j].push_back(i);
            }
        }
    }
}

std::vector<double> makeStepVector(double beg, double step, double end) {
    std::vector<double> vec;
    while (beg <= end) {
        vec.push_back(beg);
        beg += step;
    }
    return vec;
}

void removeGraphConnections(AdjacencyList &AdjacencyList, int I) {
    for (int neighbour: AdjacencyList[I]) {
        AdjacencyList[neighbour].erase(std::remove(AdjacencyList[neighbour].begin(), AdjacencyList[neighbour].end(), I),
                                       AdjacencyList[neighbour].end());
    }
    AdjacencyList[I].clear();
}

void getGraphCoordinates(const Octree &octree, const PointCloud &cloud, vector<vector<float>> &graphCoordinates,
                         AdjacencyList &AdjacencyList, int rows, int columns, float graphResolution,
                         vector<int> &realGraphIndexes) {
    pcl::PointXYZ minPt, maxPt;
    std::vector<int> pointIdxNKNSearch;
    std::vector<float> pointNKNSquaredDistance;

    pcl::getMinMax3D(*cloud, minPt, maxPt);

    auto iValues = makeStepVector(minPt.x, graphResolution, maxPt.x);
    auto jValues = makeStepVector(minPt.y, graphResolution, maxPt.y);

    float average = 0;
    for (int i = 0; i < cloud->size(); i++) {
        average += (*cloud)[i].z;
    }
    average = average / float(cloud->size());

    pcl::PointXYZ searchPoint;
    searchPoint.x = 0.0f;
    searchPoint.y = 0.0f;
    searchPoint.z = average;

    int index = 0;

    for (int j = rows - 1; j >= 0; j--) {
        for (int i = 0; i < columns; i++) {
            searchPoint.x = (float) iValues[i];
            searchPoint.y = (float) jValues[j];
            Eigen::Vector3f min_pt = Eigen::Vector3f(searchPoint.x - graphResolution * 0.4f,
                                                     searchPoint.y - graphResolution * 0.4f, -200.0f);
            Eigen::Vector3f max_pt = Eigen::Vector3f(searchPoint.x + graphResolution * 0.4f,
                                                     searchPoint.y + graphResolution * 0.4f, 200.0f);
            pcl::Indices indices;
            octree.boxSearch(min_pt, max_pt, indices);
            if (indices.empty()) {
                removeGraphConnections(AdjacencyList, index);
            } else if (indices.size() == 1) {
                graphCoordinates[index][0] = (*cloud)[indices[0]].x;
                graphCoordinates[index][1] = (*cloud)[indices[0]].y;
                graphCoordinates[index][2] = (*cloud)[indices[0]].z;
                realGraphIndexes.push_back(index);
            } else {
                float aux = (*cloud)[indices[0]].z;
                int index_aux = 0;
                for (int k = 1; k < indices.size(); k++) {
                    if ((*cloud)[indices[k]].z < aux) {
                        aux = (*cloud)[indices[k]].z;
                        index_aux = k;
                    }
                }
                graphCoordinates[index][0] = (*cloud)[indices[index_aux]].x;
                graphCoordinates[index][1] = (*cloud)[indices[index_aux]].y;
                graphCoordinates[index][2] = (*cloud)[indices[index_aux]].z;
                realGraphIndexes.push_back(index);
            }
            index++;
        }
    }
}

float distance(float p1_x, float p1_y, float p1_z, float p2_x, float p2_y, float p2_z) {
    return static_cast<float>(
            sqrt(pow(p1_x - p2_x, 2) +
            pow(p1_y - p2_y, 2) +
            pow(p1_z - p2_z, 2) * 1.0)
            );
}

float height(float p1_z, float p2_z) {
    return float(fabs(p1_z - p2_z));
}

void calculateCostList(AdjacencyList &AdjacencyList, CostList &CostList, vector<vector<float>> &graphCoordinates, int metric,
                       float distanceWeight, float heightWeight) {
    float cost;
    switch (metric) {
        case 1: // Distance
            for (int i = 0; i < AdjacencyList.size(); i++) {
                if (!AdjacencyList[i].empty()) {
                    for (int j = 0; j < AdjacencyList[i].size(); j++) {
                        cost = distance(
                                graphCoordinates[i][0],
                                graphCoordinates[i][1],
                                graphCoordinates[i][2],
                                graphCoordinates[AdjacencyList[i][j]][0],
                                graphCoordinates[AdjacencyList[i][j]][1],
                                graphCoordinates[AdjacencyList[i][j]][2]);
                        CostList[i].push_back(cost);
                    }
                }
            }
            break;
        case 2: // Height
            for (int i = 0; i < AdjacencyList.size(); i++) {
                if (!AdjacencyList[i].empty()) {
                    for (int j = 0; j < AdjacencyList[i].size(); j++) {
                        cost = height(graphCoordinates[i][2], graphCoordinates[AdjacencyList[i][j]][2]);
                        for (int k = 0; k < AdjacencyList[AdjacencyList[i][j]].size(); k++) {
                            cost += height(graphCoordinates[AdjacencyList[AdjacencyList[i][j]][k]][2],
                                           graphCoordinates[AdjacencyList[i][j]][2]);
                        }
                        CostList[i].push_back(cost);
                    }
                }
            }
            break;
        case 3: // Combined
            float dist_cost;
            float height_cost;
            for (int i = 0; i < AdjacencyList.size(); i++) {
                if (!AdjacencyList[i].empty()) {
                    for (int j = 0; j < AdjacencyList[i].size(); j++) {
                        dist_cost = distance(
                                graphCoordinates[i][0],
                                graphCoordinates[i][1],
                                graphCoordinates[i][2],
                                graphCoordinates[AdjacencyList[i][j]][0],
                                graphCoordinates[AdjacencyList[i][j]][0],
                                graphCoordinates[AdjacencyList[i][j]][2]);
                        height_cost = height(graphCoordinates[i][2], graphCoordinates[AdjacencyList[i][j]][2]);
                        for (int k = 0; k < AdjacencyList[AdjacencyList[i][j]].size(); k++) {
                            height_cost += height(graphCoordinates[AdjacencyList[AdjacencyList[i][j]][k]][2],
                                                  graphCoordinates[AdjacencyList[i][j]][2]);
                        }
                        cost = distanceWeight * dist_cost + heightWeight * height_cost;
                        CostList[i].push_back(cost);
                    }
                }
            }
            break;
        default:
            break;
    }
}

int getCorrespondentGraphNodes(vector<vector<float>> &graphCoordinates, float pointX, float pointY, vector<int> realGraphIndexes) {
    int index = 0;
    float distance = std::numeric_limits<float>::infinity();
    float temp_distance;

    for (int i = 0; i < realGraphIndexes.size(); i++) {
        temp_distance = float(
                sqrt(pow(pointX - graphCoordinates[realGraphIndexes[i]][0], 2) + pow(pointY - graphCoordinates[realGraphIndexes[i]][1], 2)));
        if (temp_distance < distance) {
            distance = temp_distance;
            index = i;
        }
    }
    return realGraphIndexes[index];
}

void getPathCoordinatesFromIndex(vector<vector<float>> &graphCoordinates, const vector<double> &pathIndexes, vector<vector<float>> &pathCoordinates) {
    for (int i = 0; i < pathIndexes.size(); i++) {
        pathCoordinates[i][0] = graphCoordinates[pathIndexes[i]][0];
        pathCoordinates[i][1] = graphCoordinates[pathIndexes[i]][1];
        pathCoordinates[i][2] = graphCoordinates[pathIndexes[i]][2];
    }
}

void printPathCoordinates(vector<vector<float>> &pathCoordinates) {
    std::cout << "Caminho gerado: " << std::endl;

    for (auto const &pathCoordinate : pathCoordinates)
        std::cout << "x: " << pathCoordinate[0] << " y: " << pathCoordinate[1] << std::endl;
}

void calculateAllCosts(const vector<double> &pathIndexes, vector<vector<float>> &graphCoordinates, std::vector<double> &allCosts,
                       float distanceWeight, float heightWeight) {
    int node1;
    int node2;
    float distanceCost = 0;
    float heightCost = 0;
    float combinedCost = 0;
    for (int i = 0; i < int(pathIndexes.size() - 1); i++) {
        node1 = int(pathIndexes[i]);
        node2 = int(pathIndexes[i + 1]);
        distanceCost += distance(
                graphCoordinates[node1][0],
                graphCoordinates[node1][1],
                graphCoordinates[node1][2],
                graphCoordinates[node2][0],
                graphCoordinates[node2][1],
                graphCoordinates[node2][2]);
        heightCost += height(graphCoordinates[node1][2], graphCoordinates[node2][2]);
        combinedCost = distanceWeight * distanceCost + heightWeight * heightCost;
    }
    allCosts[0] = distanceCost;
    allCosts[1] = heightCost;
    allCosts[2] = combinedCost;
}

int main(int argc, char **argv) {
    /*------------Setup------------*/

    string octomapFilePath {"resources/campinho10.ot"};
    int metric {1}; // 1 -> Distance, 2 -> Height, 3 -> Combined
    float distanceWeight {1.0f};
    float heightWeight {50.0f};

    auto start = high_resolution_clock::now();

    // Corredor
//    float startPositionX = 0.0f;
//    float startPositionY = 0.0f;
//    float goalPositionX = -12.10f;
//    float goalPositionY = 16.4f;

    // Campinho
    float startPositionX {13.8f};
    float startPositionY {-4.1f};
    float goalPositionX {-15.0f};
    float goalPositionY {1.25f};

    float graphAndOctreeResolution;
    PointCloud cloud = importOctomapCreatePointCloud(octomapFilePath, graphAndOctreeResolution);
    Octree octree = transformPointCloudToOctree(cloud, graphAndOctreeResolution);

    int rows {};
    int columns {};
    getMatrixDimensions(cloud, rows, columns, graphAndOctreeResolution);
    int N {rows * columns};

    vector<vector<float>> graphCoordinates(N, vector<float>(3));
    vector<double> pathIndexes;
    vector<int> realGraphIndexes;

    AdjacencyList AdjacencyList(N, std::vector<int>());
    CostList CostList(N, std::vector<float>());

    /*------------Graph Creation and Cost Calculation------------*/

    populateAdjacencyList(AdjacencyList, rows, columns);
    getGraphCoordinates(octree, cloud, graphCoordinates, AdjacencyList, rows, columns, graphAndOctreeResolution,
                        realGraphIndexes);
    calculateCostList(AdjacencyList, CostList, graphCoordinates, metric, distanceWeight, heightWeight);
    int startNodeID = getCorrespondentGraphNodes(graphCoordinates, startPositionX, startPositionY, realGraphIndexes);
    int goalNodeID = getCorrespondentGraphNodes(graphCoordinates, goalPositionX, goalPositionY, realGraphIndexes);

    std::ofstream MyFile("results/graph.txt");

    for (auto const& realGraphIndex: realGraphIndexes) {
        MyFile << graphCoordinates[realGraphIndex][0] << " "
               << graphCoordinates[realGraphIndex][1] << " "
               << graphCoordinates[realGraphIndex][2] << std::endl;
    }

    auto stop = high_resolution_clock::now();
    printExecutionTime(start, stop, "setup");

    /*------------Dijkstra------------*/

    start = high_resolution_clock::now();

    double totalCost = DijkstraPQ(AdjacencyList, CostList, startNodeID, goalNodeID, pathIndexes);

    stop = high_resolution_clock::now();
    printExecutionTime(start, stop, "dijkstra");

    cout << "Métrica utilizada: " << (metric == 1 ? "Distância" : metric == 2 ? "Altura" : "Combinada");
    cout << " - Custo total: " << totalCost << " m." << std::endl;

    /*------------Path Processing------------*/

    vector<vector<float>> pathCoordinates(pathIndexes.size(), vector<float>(3));
    getPathCoordinatesFromIndex(graphCoordinates, pathIndexes, pathCoordinates);
//    printPathCoordinates(pathCoordinates);

    /*------------Calculating costs for all metrics and saving results on file------------*/

    std::ofstream ResultsFile("results/logs.txt", std::ios_base::app);
    ResultsFile << "Octomap utilizado: " << octomapFilePath << std::endl;
    ResultsFile << "Coordenada Inicial: x = " << startPositionX << ", y = " << startPositionY << std::endl;
    ResultsFile << "Coordenada Final: x = " << goalPositionX << ", y = " << goalPositionY << std::endl;
    ResultsFile << "Métrica utilizada: " << (metric == 1 ? "Distância" : metric == 2 ? "Altura" : "Combinada")
                << std::endl;
    ResultsFile << "Peso métrica distância: " << distanceWeight << " - Peso métrica altura: " << heightWeight
                << std::endl;
    std::vector<double> allCosts(3);
    calculateAllCosts(pathIndexes, graphCoordinates, allCosts, distanceWeight, heightWeight);
    ResultsFile << "Custo total - Métrica distância: " << allCosts[0] << " m" << std::endl;
    ResultsFile << "Custo total - Métrica altura: " << allCosts[1] << " m" << std::endl;
    ResultsFile << "Custo total - Métrica combinada: " << allCosts[2] << " m" << std::endl;
    ResultsFile << "Coordenadas do caminho gerado: ";
    for (int i = 0; i < pathCoordinates.size(); i++) {
        if (i == pathCoordinates.size() - 1) {
            ResultsFile << "(" << pathCoordinates[i][0] << ", " << pathCoordinates[i][1] << ")\n";
        } else {
            ResultsFile << "(" << pathCoordinates[i][0] << ", " << pathCoordinates[i][1] << "), ";
        }
    }
    ResultsFile << "\n------------------------------------------------------\n" << std::endl;

    return 0;
}
