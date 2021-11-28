from distutils.core import setup, Extension

setup(name='cos_module',
      version='1.0',
      ext_modules=[Extension("_cos_module",
                            ["cos_module.c", "cos_module.i"])])