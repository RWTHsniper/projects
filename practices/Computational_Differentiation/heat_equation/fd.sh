#!/usr/bin/env bash


name="fd"
# Activate the licence key file by typing
export NAG_KUSARI_FILE=/home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/keys.txt

# Compile the example program in example.cpp with
g++ -I /home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/dco_cpp_v3.1.4_trial_lin64_gcc/include -c $name.cpp -o $name.o

# Link the static library part of dco/c++ (the order of object file and libary matters)
g++ -I /home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/dco_cpp_v3.1.4_trial_lin64_gcc/include $name.o /home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/dco_cpp_v3.1.4_trial_lin64_gcc/lib/libdcoc.a -o $name


