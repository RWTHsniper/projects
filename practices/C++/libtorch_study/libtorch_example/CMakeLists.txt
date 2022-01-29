# cmake -DCMAKE_PREFIX_PATH=/Users/jungjaeyong/opt/anaconda3/lib/python3.7/site-packages/torch ..

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(example-app)

set(Torch_DIR "/Users/jungjaeyong/opt/anaconda3/lib/python3.7/site-packages/torch/share/cmake/Torch/") # Torch in Python
# set(Torch_DIR "/Users/jungjaeyong/projects/practices/C++/libtorch/share/cmake/Torch/") # downloaded libtorch
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_executable(example main.cpp)
target_link_libraries(example "${TORCH_LIBRARIES}")
set_property(TARGET example PROPERTY CXX_STANDARD 17)
