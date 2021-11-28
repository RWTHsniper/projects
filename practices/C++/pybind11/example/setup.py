from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "python_example",
        sorted(glob("*.cpp")),  # Sort source files for reproducibility
        # sorted(glob("src/*.cpp")),  # Sort source files for reproducibility
    ),
]

setup(name='python_cpp_example',
    version='0.1',
    author='Benjamin Jack',
    author_email='benjamin.r.jack@gmail.com',
    description='A hybrid Python/C++ test project',
    long_description='', ext_modules=ext_modules)