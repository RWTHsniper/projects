cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")  # or 11, 14, 17, 20
project(example LANGUAGES CXX)
file(GLOB SRC
     "*.hpp"
     "*.cpp"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNO_DEBUG")

find_package(pybind11 REQUIRED)
pybind11_add_module(example ${SRC})