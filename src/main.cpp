#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/common/common.h>
#include <pcl/octree/octree_search.h>
#include <cmath>
#include <octomap/octomap.h>
#include "ros/ros.h"
#include <espeleo_control/Path.h>
#include "geometry_msgs/Polygon.h"
#include "geometry_msgs/Point32.h"
#include "geometry_msgs/PointStamped.h"
#include "std_msgs/String.h"
#include <chrono>
#include <tf/transform_listener.h>

#include "coder_array.h"
#include "Dijkstra.h"

using namespace std::chrono;

typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloud;
typedef pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> Octree;
typedef coder::array<float, 2> Matrix;
typedef std::vector<double> doubleVector;
typedef std::vector<int> intVector;
typedef std::vector<std::vector<int>> AdjacencyList;
typedef std::vector<std::vector<float>> CostList;

void
printExecutionTime(std::chrono::high_resolution_clock::time_point start,
                   std::chrono::high_resolution_clock::time_point stop,
                   const std::string &mode) {

    if (mode == "setup") {
        std::chrono::duration<double, std::milli> duration = stop - start;
        std::cout << "Setup levou "
                  << duration.count() / 1000
                  << " segundos para executar."
                  << std::endl;
    } else if (mode == "dijkstra") {
        std::chrono::duration<double, std::milli> duration = stop - start;
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

    graphAndOctreeResolution = float(octomap->getResolution());

    PointCloud cloud(new pcl::PointCloud<pcl::PointXYZ>);

    std::ofstream MyFile("/home/douglas/Área de Trabalho/Dijkstra/frontierOctomap.txt");

//    for (octomap::OcTree::leaf_bbx_iterator it = octomap->begin_leafs_bbx(octomap::point3d(10.0, -25.0, -200.0),
//                                                                          octomap::point3d(55, 50.0,
//                                                                                           200)), end = octomap->end_leafs_bbx();
    for (octomap::OcTree::leaf_iterator it = octomap->begin_leafs(), end = octomap->end_leafs();
         it != end; ++it) {
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

    rows = int(floor(fabs((double) (maxPt.y - minPt.y)) / graphResolution)) + 1;
    columns = int(floor(fabs((double) (maxPt.x - minPt.x)) / graphResolution)) + 1;
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
    vec.reserve(fabs(end - beg) / step + 1);
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

void getGraphCoordinates(const Octree &octree, const PointCloud &cloud, Matrix &graphCoordinates,
                         AdjacencyList &AdjacencyList, int rows, int columns, float graphResolution,
                         intVector &realGraphIndexes) {
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
                graphCoordinates.at(index, 0) = (*cloud)[indices[0]].x;
                graphCoordinates.at(index, 1) = (*cloud)[indices[0]].y;
                graphCoordinates.at(index, 2) = (*cloud)[indices[0]].z;
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
                graphCoordinates.at(index, 0) = (*cloud)[indices[index_aux]].x;
                graphCoordinates.at(index, 1) = (*cloud)[indices[index_aux]].y;
                graphCoordinates.at(index, 2) = (*cloud)[indices[index_aux]].z;
                realGraphIndexes.push_back(index);
            }
            index += 1;
        }
    }
//    index = 0;
////    coder::array<int, 1U> neighbours;
//    std::vector<int> neighbours;
////    coder::array<int, 1U> neighbours2;
//    std::vector<int> removeGraphNodeIndexes;
//
//    for (int j = rows - 1; j >= 0; j--) {
//        for (int i = 0; i < columns; i++) {
//            if (!AdjacencyList[index].empty()) {
//                neighbours = AdjacencyList[index];
//                for (int neighbour: neighbours) {
//                    if (fabs(graphCoordinates.at(index, 2) - graphCoordinates.at(neighbour, 2)) >
//                        graphResolution * 2.5) {
//                        removeGraphNodeIndexes.push_back(index);
//                        //                        neighbours2 = findGraphConnections(index, rows, columns);
//                        //                        for (int l = 0; l < neighbours2.size(0); l++) {
//                        //                            removeGraphNodeIndexes.push_back(neighbours2[l]);
//                        //                        }
//                    }
//                }
//            }
//            index += 1;
//        }
//    }
//
//    std::vector<int>::iterator it;
//
//    for (int removeGraphNodeIndex: removeGraphNodeIndexes) {
//        removeGraphConnections(AdjacencyList, removeGraphNodeIndex);
//        it = std::find(realGraphIndexes.begin(), realGraphIndexes.end(), removeGraphNodeIndex);
//        if (it != realGraphIndexes.end()) {
//            realGraphIndexes.erase(it);
//        }
//    }
}

float distance(float p1_x, float p1_y, float p1_z, float p2_x, float p2_y, float p2_z) {
    return float(sqrt(pow(p1_x - p2_x, 2) + pow(p1_y - p2_y, 2) + pow(p1_z - p2_z, 2) * 1.0));
}

float height(float p1_z, float p2_z) {
    return float(fabs(p1_z - p2_z));
}

void calculateCostList(AdjacencyList &AdjacencyList, CostList &CostList, Matrix &graphCoordinates, int metric,
                       float distanceWeight, float heightWeight) {
    float cost;
    switch (metric) {
        case 1: // Distance
            for (int i = 0; i < AdjacencyList.size(); i++) {
                if (!AdjacencyList[i].empty()) {
                    for (int j = 0; j < AdjacencyList[i].size(); j++) {
                        cost = distance(
                                graphCoordinates.at(i, 0),
                                graphCoordinates.at(i, 1),
                                graphCoordinates.at(i, 2),
                                graphCoordinates.at(AdjacencyList[i][j], 0),
                                graphCoordinates.at(AdjacencyList[i][j], 1),
                                graphCoordinates.at(AdjacencyList[i][j], 2));
                        CostList[i].push_back(cost);
                    }
                }
            }
            break;
        case 2: // Height
            for (int i = 0; i < AdjacencyList.size(); i++) {
                if (!AdjacencyList[i].empty()) {
                    for (int j = 0; j < AdjacencyList[i].size(); j++) {
                        cost = height(graphCoordinates.at(i, 2), graphCoordinates.at(AdjacencyList[i][j], 2));
                        for (int k = 0; k < AdjacencyList[AdjacencyList[i][j]].size(); k++) {
                            cost += height(graphCoordinates.at(AdjacencyList[AdjacencyList[i][j]][k], 2),
                                           graphCoordinates.at(AdjacencyList[i][j], 2));
                        }
//                        for (int k = 0; k < AdjacencyList[AdjacencyList[i][j]].size(); k++) {
//                            cost += height(graphCoordinates.at(AdjacencyList[AdjacencyList[i][j]][k], 2),
//                                           graphCoordinates.at(i, 2));
//                        }
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
                                graphCoordinates.at(i, 0),
                                graphCoordinates.at(i, 1),
                                graphCoordinates.at(i, 2),
                                graphCoordinates.at(AdjacencyList[i][j], 0),
                                graphCoordinates.at(AdjacencyList[i][j], 1),
                                graphCoordinates.at(AdjacencyList[i][j], 2));
                        height_cost = height(graphCoordinates.at(i, 2), graphCoordinates.at(AdjacencyList[i][j], 2));
                        for (int k = 0; k < AdjacencyList[AdjacencyList[i][j]].size(); k++) {
                            height_cost += height(graphCoordinates.at(AdjacencyList[AdjacencyList[i][j]][k], 2),
                                                  graphCoordinates.at(AdjacencyList[i][j], 2));
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

int getCorrespondentGraphNodes(Matrix &graphCoordinates, float pointX, float pointY, intVector realGraphIndexes) {
    int index = 0;
    float distance = std::numeric_limits<float>::infinity();
    float temp_distance;

    for (int i = 0; i < realGraphIndexes.size(); i++) {
        temp_distance = float(
                sqrt(pow(pointX - graphCoordinates.at(realGraphIndexes[i], 0), 2) +
                     pow(pointY - graphCoordinates.at(realGraphIndexes[i], 1), 2)));
        if (temp_distance < distance) {
            distance = temp_distance;
            index = i;
        }
    }
    return realGraphIndexes[index];
}

void getPathCoordinatesFromIndex(Matrix &graphCoordinates, const doubleVector &pathIndexes, Matrix &pathCoordinates) {

    for (int i = 0; i < int(pathIndexes.size()); i++) {
        pathCoordinates.at(i, 0) = graphCoordinates.at(int(pathIndexes[i]), 0);
        pathCoordinates.at(i, 1) = graphCoordinates.at(int(pathIndexes[i]), 1);
        pathCoordinates.at(i, 2) = graphCoordinates.at(int(pathIndexes[i]), 2);
    }
}

void printPathCoordinates(Matrix &pathCoordinates) {
    std::cout << "Caminho gerado: " << std::endl;

    for (int i = 0; i < pathCoordinates.size(0); i++) {
        std::cout << "x: " << pathCoordinates.at(i, 0) << " y: " << pathCoordinates.at(i, 1) << std::endl;
    }
}

geometry_msgs::Polygon createPathMessageRos(Matrix &pathCoordinates) {
    geometry_msgs::Polygon path_msg;
    geometry_msgs::Point32 p;
    p.z = 0.0;

    std::ofstream MyFile("/home/douglas/Área de Trabalho/Dijkstra/path.txt");

    for (int i = 0; i < pathCoordinates.size(0); i++) {
        p.x = pathCoordinates.at(i, 0);
        p.y = pathCoordinates.at(i, 1);
        path_msg.points.push_back(p);

        MyFile << pathCoordinates.at(i, 0) << " "
               << pathCoordinates.at(i, 1) << " "
               << pathCoordinates.at(i, 2) << std::endl;
    }

    return path_msg;
}

void calculateAllCosts(const doubleVector &pathIndexes, Matrix &graphCoordinates, std::vector<double> &allCosts,
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
                graphCoordinates.at(node1, 0),
                graphCoordinates.at(node1, 1),
                graphCoordinates.at(node1, 2),
                graphCoordinates.at(node2, 0),
                graphCoordinates.at(node2, 1),
                graphCoordinates.at(node2, 2));
        heightCost += height(graphCoordinates.at(node1, 2), graphCoordinates.at(node2, 2));
        combinedCost = distanceWeight * distanceCost + heightWeight * heightCost;
    }
    allCosts[0] = distanceCost;
    allCosts[1] = heightCost;
    allCosts[2] = combinedCost;
}

int main(int argc, char **argv) {
    /*------------ROS Setup------------*/

    ros::init(argc, argv, "dijkstra_planning");
    ros::NodeHandle n;
    ros::Duration(0.5).sleep();

    std::string octomapFilePath; // Arquivo octomap a ser importado
    int metric; // Métrica utilizada no planejamento de caminhos
    float distanceWeight;
    float heightWeight;
    std::string sourceFrame; // Frames para pegar pose do robô
    std::string targetFrame;

    try {
        n.param("/dijkstra_planning/octomapFilePath", octomapFilePath, std::string("tmp/"));
        n.param("/dijkstra_planning/metric", metric, int(1));
        n.param("/dijkstra_planning/distanceWeight", distanceWeight, float(1.0));
        n.param("/dijkstra_planning/heightWeight", heightWeight, float(50.0));
        n.param("/dijkstra_planning/sourceFrame", sourceFrame, std::string("base_init"));
        n.param("/dijkstra_planning/targetFrame", targetFrame, std::string("ground_truth"));
    } catch (int e) {
        ROS_INFO(
                "\033[1;31m---->\033[0m Exception occurred when importing parameters in Dijkstra Planning node. Exception Nr. %d",
                e);
    }

    // Tópico onde será publicado o caminho gerado
    ros::Publisher pub_traj = n.advertise<espeleo_control::Path>("espeleo/traj_points", 1);
    espeleo_control::Path path;
    path.closed_path_flag = false;
    path.insert_n_points = 3;
    path.filter_path_n_average = 3;


    // Obtendo coordenada inicial
//    tf::TransformListener listener;
//    tf::StampedTransform transform;
//    std::cout << "Esperando coordenada inicial do caminho ..." << std::endl;
//    int tf_read = 0;
//    while (tf_read == 0) {
//        try {
//            listener.lookupTransform(sourceFrame, targetFrame, ros::Time(0), transform);
//            tf_read = 1;
//        } catch (tf::TransformException &ex) {
//            ros::Duration(1.0).sleep();
//        }
//    }
//    std::cout << "Coordenada inicial recebida -> x: " << float(transform.getOrigin().x()) << " y: "
//              << float(transform.getOrigin().y()) << std::endl;

//    std::cout << "Esperando coordenada inicial do caminho ..." << std::endl;
//    geometry_msgs::PointStampedConstPtr startPosePointer = ros::topic::waitForMessage<geometry_msgs::PointStamped>(
//            "clicked_point");
//    geometry_msgs::PointStamped startPose = *startPosePointer;
//    std::cout << "Coordenada final recebida -> x: " << float(startPose.point.x) << " y: " << float(startPose.point.y)
//              << std::endl;


    // Obtendo coordenada final
//    std::cout << "Esperando coordenada final do caminho ..." << std::endl;
//    geometry_msgs::PointStampedConstPtr goalPosePointer = ros::topic::waitForMessage<geometry_msgs::PointStamped>(
//            "clicked_point");
//    geometry_msgs::PointStamped goalPose = *goalPosePointer;
//    std::cout << "Coordenada final recebida -> x: " << float(goalPose.point.x) << " y: " << float(goalPose.point.y)
//              << std::endl;


    /*------------Setup------------*/
    auto start = high_resolution_clock::now();

//    auto startPositionX = float(transform.getOrigin().x());
//    auto startPositionY = float(transform.getOrigin().y());
//    auto goalPositionX = float(goalPose.point.x);
//    auto goalPositionY = float(goalPose.point.y);
//
//    auto startPositionX = float(startPose.point.x);
//    auto startPositionY = float(startPose.point.y);

// Campinho
//    float startPositionX = 13.8f;
//    float startPositionY = -4.1f;
//    float goalPositionX = -15.0f;
//    float goalPositionY = 1.25f;

// Corredor
//    float startPositionX = 0.0f;
//    float startPositionY = 0.0f;
//    float goalPositionX = -12.10f;
//    float goalPositionY = 16.4f;

// Corredor UFMG 1
//    float startPositionX = 17.1577f;
//    float startPositionY = -0.292194f;
//    float goalPositionX = -0.477881f;
//    float goalPositionY = -6.2912f;

// Corredor UFMG 2
//    float startPositionX = -0.597627f;
//    float startPositionY = -6.46333f;
//    float goalPositionX = -8.75299f;
//    float goalPositionY = -0.400397f;

// Corredor UFMG 3
//    float startPositionX = -8.83121f;
//    float startPositionY = -0.641502f;
//    float goalPositionX = 1.66436f;
//    float goalPositionY = 3.62138f;

// Corredor UFMG 4
    float startPositionX = 1.72808f;
    float startPositionY = 3.48767f;
    float goalPositionX = 14.5491f;
    float goalPositionY = -0.429702f;

    float graphAndOctreeResolution; // Resolução da octree e grafo em metros

    PointCloud cloud = importOctomapCreatePointCloud(octomapFilePath, graphAndOctreeResolution);
    Octree octree = transformPointCloudToOctree(cloud, graphAndOctreeResolution);

    int rows = 0;
    int columns = 0;
    getMatrixDimensions(cloud, rows, columns, graphAndOctreeResolution);
    int N = rows * columns;

    Matrix graphCoordinates;
    Matrix pathCoordinates;
    doubleVector pathIndexes;
    intVector realGraphIndexes;

    AdjacencyList AdjacencyList(N, std::vector<int>(0));
    CostList CostList(N, std::vector<float>(0));

    graphCoordinates.set_size(N, 3);

    /*------------Graph Creation and Cost Calculation------------*/

    populateAdjacencyList(AdjacencyList, rows, columns);
    getGraphCoordinates(octree, cloud, graphCoordinates, AdjacencyList, rows, columns, graphAndOctreeResolution,
                        realGraphIndexes);
    calculateCostList(AdjacencyList, CostList, graphCoordinates, metric, distanceWeight, heightWeight);
    int startNodeID = getCorrespondentGraphNodes(graphCoordinates, startPositionX, startPositionY, realGraphIndexes);
    int goalNodeID = getCorrespondentGraphNodes(graphCoordinates, goalPositionX, goalPositionY, realGraphIndexes);

    std::ofstream MyFile("/home/douglas/Área de Trabalho/Dijkstra/graph.txt");

    for (int realGraphIndex: realGraphIndexes) {
        MyFile << graphCoordinates.at(realGraphIndex, 0) << " "
               << graphCoordinates.at(realGraphIndex, 1) << " "
               << graphCoordinates.at(realGraphIndex, 2) << std::endl;
    }

    auto stop = high_resolution_clock::now();
    printExecutionTime(start, stop, "setup");

    /*------------Dijkstra------------*/

    start = high_resolution_clock::now();

    double totalCost = Dijkstra(AdjacencyList, CostList, startNodeID, goalNodeID, pathIndexes);

    stop = high_resolution_clock::now();
    printExecutionTime(start, stop, "dijkstra");

    std::cout << "Métrica utilizada: " << (metric == 1 ? "Distância" : metric == 2 ? "Altura" : "Combinada");
    std::cout << " - Custo total: " << totalCost << " m." << std::endl;

    /*------------Path Processing------------*/

    pathCoordinates.set_size(int(pathIndexes.size()), 3);
    getPathCoordinatesFromIndex(graphCoordinates, pathIndexes, pathCoordinates);
//    printPathCoordinates(pathCoordinates);
    path.path = createPathMessageRos(pathCoordinates);

    /*------------Calculating costs for all metrics and saving results on file------------*/

    std::ofstream ResultsFile("/home/douglas/catkin_ws/src/dijkstra_planning/resources/logs.txt", std::ios_base::app);
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
    for (int i = 0; i < pathCoordinates.size(0); i++) {
        if (i == pathCoordinates.size(0) - 1) {
            ResultsFile << "(" << pathCoordinates.at(i, 0) << ", " << pathCoordinates.at(i, 1) << ")\n";
        } else {
            ResultsFile << "(" << pathCoordinates.at(i, 0) << ", " << pathCoordinates.at(i, 1) << "), ";
        }
    }
    ResultsFile << "\n------------------------------------------------------\n" << std::endl;

    /*------------Path Publishing------------*/

    while (ros::ok()) {
        if (pub_traj.getNumSubscribers() > 0) {
            ros::Duration(0.5).sleep();
            std::cout << "Enviando caminho gerado para controle do robô.\n" << std::endl;
            pub_traj.publish(path);
            break;
        }
    }

    return 0;
}
