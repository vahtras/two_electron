"""Some docstring"""
import os
import pathlib

import numpy as np
import two.eri


class TestTwo:
    """Two-fock test suite"""

    def setup_method(self):
        """Setup supporting directory"""
        name, _ = os.path.splitext(__file__)
        self.suppdir = pathlib.Path(name + ".d")
        self.aotwoint = self.suppdir/"AOTWOINT"
        self.reader = two.eri.Reader(self.aotwoint)

        self.faref = np.array([
            [2.02818057, 0.26542036, 0.00000000, 0.06037429, 0.00000000, 0.00000000],
            [0.26542036, 0.74551226, 0.00000000, 0.28646605, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 1.01061061, -0.47343699, 0.00000000, 0.00000000],
            [0.06037429, 0.28646605, -0.47343699, 0.86039561, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.01061061, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.01061061],
            ])

        self.fbref = np.array([
            [2.07812203, 0.35828051, 0.00000000, 0.09571548, 0.00000000, 0.00000000],
            [0.35828051, 1.03607456, 0.00000000, 0.40151695, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 1.07479494, -0.50700370, 0.00000000, 0.00000000],
            [0.09571548, 0.40151695, -0.50700370, 0.93374109, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.07479494, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.07479494],
            ])

    def test_fab_p(self):
        """Test alpha and beta Fock matrix"""
        d_a = np.zeros((6, 6))
        d_b = np.zeros((6, 6))
        d_a[0, 0] = 1
        d_a[1, 1] = 1
        d_b[0, 0] = 1
        (f_a, f_b), = self.reader.fockab((d_a, d_b), f2py=False)
        np.testing.assert_allclose(f_a, self.faref)
        np.testing.assert_allclose(f_b, self.fbref)

    def test_fab_f(self):
        """Test alpha and beta Fock matrix, Fortran version"""
        d_a = np.zeros((6, 6))
        d_b = np.zeros((6, 6))
        d_a[0, 0] = 1
        d_a[1, 1] = 1
        d_b[0, 0] = 1
        (f_a, f_b), = self.reader.fockab((d_a, d_b), f2py=True)
        np.testing.assert_allclose(f_a, self.faref)
        np.testing.assert_allclose(f_b, self.fbref)

    def test_f_p(self):
        "Test total Fock, Python version"""
        d_a = np.zeros((6, 6))
        d_b = np.zeros((6, 6))
        d_a[0, 0] = 1
        d_a[1, 1] = 1
        d_b[0, 0] = 1
        dtot = d_a + d_b
        ftot = self.reader.fock(dtot, f2py=False)

        fref = 0.5*(self.faref+self.fbref)
        np.testing.assert_allclose(ftot, fref)

    def test_f_f(self):
        "Test total Fock, Fortran version"""
        d_a = np.zeros((6, 6))
        d_b = np.zeros((6, 6))
        d_a[0, 0] = 1
        d_a[1, 1] = 1
        d_b[0, 0] = 1
        dtot = d_a + d_b
        ftot = self.reader.fock(dtot, f2py=True)

        fref = 0.5*(self.faref+self.fbref)
        np.testing.assert_allclose(ftot, fref)

    def test_fs_p(self):
        "Test spin Fock, Python version"""
        d_a = np.zeros((6, 6))
        d_b = np.zeros((6, 6))
        d_a[0, 0] = 1
        d_a[1, 1] = 1
        d_b[0, 0] = 1
        dspin = d_a - d_b
        fspin = self.reader.fock(dspin, hfc=0, f2py=False)

        fref = 0.5*(self.faref-self.fbref)

        np.testing.assert_allclose(fspin, fref, atol=1e-8)

    def test_fs_f(self):
        "Test spin Fock, Python version"""
        d_a = np.zeros((6, 6))
        d_b = np.zeros((6, 6))
        d_a[0, 0] = 1
        d_a[1, 1] = 1
        d_b[0, 0] = 1
        dspin = d_a - d_b
        fspin = self.reader.fock(dspin, hfc=0, f2py=True)

        fref = 0.5*(self.faref-self.fbref)

        np.testing.assert_allclose(fspin, fref, atol=1e-8)
