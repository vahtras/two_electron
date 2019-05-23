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
    name='two_electron',
    packages=['two'],
    install_requires=["blocked-matrix-utils", "fortran-binary"],
    ext_modules=[ext]
)
