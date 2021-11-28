from distutils.core import setup, Extension
import numpy

setup(name="cos_doubles",
ext_modules=[Extension("_cos_doubles",
      sources=["cos_doubles.c", "cos_doubles.i"],
      include_dirs=[numpy.get_include()])])