from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# Explicit Extension objects so module names are 'calc', 'parse', 'thread_wrapper'
# (not 'src.calc' etc.), keeping afcn.py's bare imports working.
ext_modules = cythonize(
    [
        Extension("calc", sources=["src/calc.pyx"]),
        Extension("parse", sources=["src/parse.pyx"]),
        Extension("thread_wrapper", sources=["src/thread_wrapper.pyx"]),
    ],
    compiler_directives={"language_level": "3"},
)

setup(
    ext_modules=ext_modules,
    include_dirs=[numpy.get_include()],
    scripts=["src/afcn.py", "src/convert_genexpc_to_afcn.py"],
)
