cmake_minimum_required(VERSION 3.0.0)
project(
  MultiAsset
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

include(CTest)
enable_testing()

add_executable(${PROJECT_NAME} main.cpp utils.cpp sde.cpp MyTensor.cpp)
# set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -fopenmp -DNO_DEBUG")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/eigen/) # Source codes are there

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)
