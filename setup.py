# setuptools is not used here but picked up by numpy.distutils.core.setup
# if setuptools is in sys.modules final call will be to setuptools.setup

import setuptools
import sys


try:
    from numpy.distutils.core import setup, Extension
except ModuleNotFoundError:
    import subprocess
    subprocess.call(f"{sys.executable} -m pip install numpy".split())
    from numpy.distutils.core import setup, Extension

__version__ = "1.1.1"

ext = Extension(
    name='sirfck',
    sources=['two/sirfck.f'],
    libraries=['blas']
)

setup(
    name='two_electron',
    version=__version__,
    author="Olav Vahtras",
    author_email="olav.vahtras@gmail.com",
    url="https://github.com/vahtras/two_electron",
    install_requires=["blocked-matrix-utils", "fortran-binary"],
    packages=['two'],
    ext_modules=[ext]
)
