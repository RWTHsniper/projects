cmake_minimum_required(VERSION 3.0.0)
# project(g++ VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(MultiAsset main.cpp utils.cpp sde.cpp MyTensor.cpp)
# set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -fopenmp -DNO_DEBUG")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
