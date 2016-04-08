try:
    import numpy
except ImportError:
    import subprocess
    subprocess.call("pip install numpy", shell=True)

from numpy.distutils.core import setup, Extension
ext = Extension(
    name='sirfck',
    sources=['sirfck.f'],
    libraries=['blas']
)

setup(
    name='two',
    ext_modules=[ext]
)
