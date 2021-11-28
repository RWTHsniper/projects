/*
c++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)

$ python
Python 2.7.10 (default, Aug 22 2015, 20:33:39)
[GCC 4.2.1 Compatible Apple LLVM 7.0.0 (clang-700.0.59.1)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
import example
example.add(1, 2)
3L
*/

#include "example.hpp"

namespace py = pybind11;


int add(int i, int j) {
    return i + j;
}

struct Pet {
    Pet(const std::string &name) : name(name) { }
    void setName(const std::string &name_) { name = name_; }
    const std::string &getName() const { return name; }

    std::string name;
};

py::array_t<double> add_arrays_test(py::array_t<double> input1, py::array_t<double> input2) { // py::buffer handles general python object, but py::array_t only serves numpy array

    auto start = std::chrono::steady_clock::now();
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();
    
    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /* No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

    auto end = std::chrono::steady_clock::now();
 
    std::cout << "Elapsed time for initialization in nanoseconds: "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << " ns" << std::endl;

    for (size_t idx = 0; idx < buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];
    // add_arrays(ptr3, ptr1, ptr2, static_cast<size_t>(buf1.shape[0])); // My reusable C++ function

    end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time for completion in nanoseconds: "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << " ns" << std::endl;

    return result;
}


py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2) { // py::buffer handles general python object, but py::array_t only serves numpy array
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();

    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /* No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

    for (size_t idx = 0; idx < buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];
    // add_arrays(ptr3, ptr1, ptr2, static_cast<size_t>(buf1.shape[0])); // My reusable C++ function

    return result;
}

void add_arrays(py::array_t<double> result, py::array_t<double> input1, py::array_t<double> input2) {
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();
    py::buffer_info buf3 = result.request();

    if (buf1.ndim != 1 || buf2.ndim != 1 || buf3.ndim != 1 )
        throw std::runtime_error("Number of dimensions must be one");

    if (!(buf1.size == buf2.size && buf2.size == buf3.size && buf3.size == buf1.size)) {
        throw std::runtime_error("Input shapes must match");
    }

    /* No pointer is passed, so NumPy will allocate the buffer */

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

    for (size_t idx = 0; idx < buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];
    // add_arrays(ptr3, ptr1, ptr2, static_cast<size_t>(buf1.shape[0])); // My reusable C++ function

}

void add_arrays(double *out, double *ptr1, double *ptr2, size_t n){
    for (size_t idx = 0; idx < n; idx++)
    out[idx] = ptr1[idx] + ptr2[idx];
}

void add_arrays_D(py::array_t<double> result, py::array_t<double> input1, py::array_t<double> input2) {
    auto i1 = result.unchecked<1>(); // Will throw if ndim != 1 or flags.writeable is false
    auto i2 = result.unchecked<1>(); // Will throw if ndim != 1 or flags.writeable is false
    auto r = result.mutable_unchecked<1>(); // Will throw if ndim != 1 or flags.writeable is false
    for (py::ssize_t i = 0; i < r.shape(0); i++)
        r(i) = i1(i) + i2(i);
}


PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    // regular notation
    m.def("add1", &add, py::arg("i"), py::arg("j"));
    m.def("add", &add, "A function which adds two numbers",
      py::arg("i"), py::arg("j"));
    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string &>())
        .def("setName", &Pet::setName)
        .def_readwrite("name", &Pet::name) // directly expose name to Python. It can be directly assess from Python.
        .def("__repr__",
            [](const Pet &a) {
                return "<example.Pet named '" + a.name + "'>";
            }
        ) // print(p)
        .def("getName", &Pet::getName);
    // m.def("add_arrays", &add_arrays, "Add two NumPy arrays"); // skipped describing argument types
    m.def("add_arrays", py::overload_cast<py::array_t<double>, py::array_t<double>>(&add_arrays)); // specified argument types because of overloading
    m.def("add_arrays_test", &add_arrays_test);
    m.def("add_arrays", py::overload_cast<py::array_t<double>, py::array_t<double>, py::array_t<double>>(&add_arrays));
    m.def("add_arrays_D", py::overload_cast<py::array_t<double>, py::array_t<double>, py::array_t<double>>(&add_arrays));

}


/*
Taking keyword arguments
m.def("add", &add, "A function which adds two numbers",
      py::arg("i"), py::arg("j"));

import example
example.add(i=1, j=2)
3L

% python
import example
p = example.Pet("Molly")
print(p)
<example.Pet object at 0x10cd98060>
p.getName()
u'Molly'
p.setName("Charly")
p.getName()
u'Charly'

*/