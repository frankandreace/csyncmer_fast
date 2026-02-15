import shutil
from pathlib import Path
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ROOT_DIR = Path(__file__).parent.resolve()
REPO_ROOT = ROOT_DIR.parent.parent
HEADER = 'csyncmer_fast.h'
LOCAL_HEADER = ROOT_DIR / HEADER

# Copy header from repo root into package dir so it's available during
# sdist/wheel builds (which run from a temp directory).
if not LOCAL_HEADER.exists():
    repo_header = REPO_ROOT / HEADER
    if repo_header.exists():
        shutil.copy2(repo_header, LOCAL_HEADER)

if not LOCAL_HEADER.exists():
    raise RuntimeError(
        f"Cannot find {HEADER}. Build from the repository or ensure "
        f"{HEADER} is present in {ROOT_DIR}."
    )

ext_modules = [
    Pybind11Extension(
        "csyncmer_fast._bindings",
        sources=[
            "csyncmer_fast/_bindings.cpp",
        ],
        include_dirs=[
            str(ROOT_DIR),
        ],
        cxx_std=17,
        extra_compile_args=[
            '-O3',
            '-march=native',
            '-shared',
            '-fPIC',
            '-ffast-math',
            '-funroll-loops',
            '-flto',
            '-fprefetch-loop-arrays',
            '-mavx2',
        ],
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
