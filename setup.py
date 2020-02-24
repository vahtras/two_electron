import setuptools

try:
    import numpy
except ImportError:
    import subprocess
    subprocess.call("pip install numpy", shell=True)

from numpy.distutils.core import setup, Extension
ext = Extension(
    name='sirfck',
    sources=['two/sirfck.f'],
    libraries=['blas']
)

setup(
    author="Olav Vahtras",
    author_email="olav.vahtras@gmail.com",
    version='1.0.0',
    url='https://github.com/vahtras/two_electron',
    name='two-electron',
    packages=['two'],
    install_requires=["blocked-matrix-utils", "fortran-binary"],
    ext_modules=[ext]
)
