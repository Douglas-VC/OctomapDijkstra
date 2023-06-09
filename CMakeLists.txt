cmake_minimum_required(VERSION 3.5.1)
project(dijkstra_planning)
set(CMAKE_CXX_STANDARD 11)

find_package(PCL 1.7 REQUIRED COMPONENTS common io)
find_package(octomap REQUIRED)

include_directories(
        include
        ${PCL_INCLUDE_DIRS}
        ${OCTOMAP_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

SET(TARGET_H src/Dijkstra.h)

SET(TARGET_SRC src/Dijkstra.cpp)

add_executable(DijkstraPathPlanning src/main.cpp ${TARGET_H} ${TARGET_SRC})

target_link_libraries(DijkstraPathPlanning ${PCL_LIBRARIES} ${OCTOMAP_LIBRARIES})