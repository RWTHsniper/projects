#include <stdlib.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <chrono>

void add_arrays(double *out, double *ptr1, double *ptr2, size_t n);
