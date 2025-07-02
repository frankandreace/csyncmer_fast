import os
import sys
from pathlib import Path
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


def check_nthash_library(lib_dirs):
    for lib_dir in lib_dirs:
        potential_path1 = os.path.join(lib_dir, 'libnthash.a') 
        potential_path2 = os.path.join(lib_dir, 'lib/libnthash.a') 
        if os.path.exists(potential_path1) or os.path.exists(potential_path2):
            return True
    return False

library_dirs=['/usr/local/nthash','/usr/local/lib/nthash','/usr/lib/nthash','/usr/local/lib','/usr/lib/']

# CHECK THAT NTHASH LIBRARY IS PRESENT
ROOT_DIR = Path(__file__).parent.resolve()
print(f">>>ROOT DIR IS: {ROOT_DIR}")

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
    Pybind11Extension(
        "csyncmer_fast._bindings",
        sources=[
            "src/csyncmer_fast/_bindings.cpp",
            "src/csyncmer_fast/c_lib/vendor/syng/utils_d.c",
            "src/csyncmer_fast/c_lib/vendor/syng/seqhash.c",
        ],
        include_dirs=[
            'src/csyncmer_fast/c_lib/include',
            'src/csyncmer_fast/c_lib/vendor',
        ],
        cxx_std=17,
        extra_compile_args = [
            '-O3',
            '-march=native',
            '-shared',
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
    cmdclass={"build_ext":build_ext},
    # Make sure to include package data
    package_data={
        "csyncmer_fast": [
            "c_lib/dev/*.h",
            "c_lib/vendor/syng/*.h",
            "c_lib/vendor/nthash.hpp",
            "c_lib/include/csyncmer_fast/*.h",
        ],
    },
)