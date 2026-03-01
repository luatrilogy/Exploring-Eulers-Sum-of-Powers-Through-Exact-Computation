from setuptools import setup, Extension
import pybind11

boost_include = r"A:\My Stuff\Projects\Python Stuff\Euler Sum of Powers Conjecture\vcpkg\installed\x64-windows\include"

ext_modules = [
    Extension(
        "euler_ext",
        ["euler_ext.cpp"],
        include_dirs=[pybind11.get_include(), boost_include],
        language="c++",
        extra_compile_args=["/O2"],
    )
]

setup(
    name="euler_ext",
    version="0.0.1",
    ext_modules=ext_modules,
)