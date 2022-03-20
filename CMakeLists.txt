cmake_minimum_required(VERSION 3.5)

project(quicksort LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_CXX_FLAGS "-fsanitize=thread")

find_package(TBB REQUIRED)

add_executable(quicksort main.cpp)
add_executable(bfs bfs.cpp)

target_link_libraries(quicksort PRIVATE TBB::tbb)
target_link_libraries(bfs PRIVATE TBB::tbb)
