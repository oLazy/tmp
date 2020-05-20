cmake_minimum_required(VERSION 3.16)
project(mt1danisModelParams)
find_package(Boost REQUIRED serialization)
set(CMAKE_CXX_STANDARD 14)

add_executable(mt1danisModelParams main.cpp)
target_link_libraries(mt1danisModelParams ${Boost_LIBRARIES})