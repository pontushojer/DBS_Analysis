from distutils.core import setup
from Cython.Build import cythonize
setup(ext_modules=cythonize(["dbs_analysis/hamming_cython_solution.pyx"]))

# adapted from https://github.com/kwmsmith/scipy-2015-cython-tutorial
