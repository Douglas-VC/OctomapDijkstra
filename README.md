# OctomapDijkstra

Path Planning implementation using the Octomap library and Dijkstra's algorithm. The program needs a octomap(.ot or .bt file, built from point cloud data, for example), which is processed and then used to create a graph. Dijsktra's algorithm is then applied to find the minimum cost path between start and end points, returning total cost and the path generated. Includes multiple implementations of the Dijkstra's algorithm, one using a priority_queue and the other using a set.
