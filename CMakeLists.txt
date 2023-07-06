cmake_minimum_required(VERSION 3.24)
project(OctomapDijkstra)
set(CMAKE_CXX_STANDARD 17)

set(octomap_DIR "/usr/local/octomap-1.9.6/octomap/lib/cmake/octomap/")

find_package(PCL 1.12 REQUIRED COMPONENTS common io)
find_package(octomap REQUIRED)

include_directories(include ${PCL_INCLUDE_DIRS} ${OCTOMAP_INCLUDE_DIRS})

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

SET(TARGET_H src/Dijkstra.h)
SET(TARGET_SRC src/Dijkstra.cpp src/DijkstraPQ.cpp src/DijkstraPQ.h)

add_executable(OctomapDijkstra src/main.cpp ${TARGET_H} ${TARGET_SRC})
target_link_libraries(OctomapDijkstra ${PCL_LIBRARIES} ${OCTOMAP_LIBRARIES})