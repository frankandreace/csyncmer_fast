import os
import sys
from setuptools import setup, Extension
import pybind11

def check_nthash_library(lib_dirs):
    for lib_dir in lib_dirs:
        potential_path1 = os.path.join(lib_dir, 'libnthash.a') 
        potential_path2 = os.path.join(lib_dir, 'lib/libnthash.a') 
        if os.path.exists(potential_path1) or os.path.exists(potential_path2):
            return True
    return False

library_dirs=['/usr/local/nthash','/usr/local/lib/nthash','/usr/lib/nthash','/usr/local/lib','/usr/lib/']

print(os.path.abspath('../include'))
# CHECK THAT NTHASH LIBRARY IS PRESENT
if not check_nthash_library(library_dirs):
    sys.stderr.write(
        "Error: nthash library not found in specified directories.\n"
        "Please ensure nthash is installed and provide its path in library_dirs.\n"
        "Library dirs are:"
    )
    for library_dir in library_dirs:
        sys.stderr.write(f" {library_dir}")
    sys.stderr.write("\n")
    sys.exit(1)

ext_modules = [
    Extension('csyncmer_fast._bindings',
              sources=['csyncmer_fast/_bindings.cpp'],
              include_dirs=[
                  pybind11.get_include(),
                  '../include'  # Changed path
              ],
              language='c++',
              extra_compile_args =  ['-O3',
                                    '-march=native',
                                    '-shared',
                                    '-std=c++17',
                                    '-fPIC',
                                    '-ffast-math',
                                    '-funroll-loops',
                                    '-flto',
                                    '-fprefetch-loop-arrays',
                                    '-mavx2',
                                    '-v',
                                    ],
              libraries=['nthash'],
              library_dirs=library_dirs,
              ),
]


setup(
    ext_modules=ext_modules,
)