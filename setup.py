from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.1"

INSTALL_REQUIRES = [
    'rdkit~=2023.9.1',
    'numpy~=1.26.3',
    'pybind11~=2.11.1'
]
PACKAGES = [
    'gesim'
]

ext_modules = [
    Pybind11Extension(
        "gesim.convert",
        ["src/convert.cpp"],
        extra_compile_args=['-O3', '-Wall', '-std=c++11'],
        define_macros=[("VERSION_INFO", __version__)],
    ),
    Pybind11Extension(
        "gesim.gdb",
        ["src/gdb.cpp"],
        extra_compile_args=['-O3', '-Wall', '-std=c++11'],
        define_macros=[("VERSION_INFO", __version__)],
    ),
    Pybind11Extension(
        "gesim.graph_entropy",
        ["src/graph_entropy.cpp"],
        extra_compile_args=['-O3', '-Wall', '-std=c++11'],
        define_macros=[("VERSION_INFO", __version__)],
    ),
]

setup(
    name="gesim",
    version=__version__,
    author="Hiroaki Shiokawa",
    author_email="",
    url="",
    description="graph entropy similarity",
    long_description="",
    ext_modules=ext_modules,
    #extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.11",
    install_requires=INSTALL_REQUIRES,
    packages=PACKAGES,
)