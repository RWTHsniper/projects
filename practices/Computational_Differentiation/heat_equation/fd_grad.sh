#!/usr/bin/env bash


name="fd_grad"

export NAG_KUSARI_FILE=/home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/keys.txt

g++ -I /home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/dco_cpp_v3.1.4_trial_lin64_gcc/include -c $name.cpp -o $name.o

g++ -I /home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/dco_cpp_v3.1.4_trial_lin64_gcc/include $name.o /home/jaeyong/Desktop/studies/16ws/Computational_Differentiation/dco_cpp_v3.1.4_trial_lin64_gcc/lib/libdcoc.a -o $name


