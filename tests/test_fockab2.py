import os
import pathlib

import numpy as np
from two import fockab, fock


class TestTwo:
    """Two-fock test suite"""

    def setup_method(self):
        """Setup supporting directory"""
        self.suppdir = pathlib.Path(__file__).with_suffix(".d")
        self.aotwoint = self.suppdir/"AOTWOINT"

        self.faref = np.loadtxt(self.suppdir/'fa.ref')
        self.fbref = np.loadtxt(self.suppdir/'fb.ref')

        self.daref = np.loadtxt(self.suppdir/'da.ref')
        self.dbref = np.loadtxt(self.suppdir/'db.ref')

    def test_fab_p(self):
        """Test alpha and beta Fock matrix"""
        d_a, d_b = self.daref, self.dbref
        (f_a, f_b), = fockab((d_a, d_b), filename=self.aotwoint, f2py=False)
        np.testing.assert_allclose(f_a, self.faref)
        np.testing.assert_allclose(f_b, self.fbref)

    def test_fab_f(self):
        """Test alpha and beta Fock matrix, Fortran version"""
        d_a, d_b = self.daref, self.dbref
        (f_a, f_b), = fockab((d_a, d_b), filename=self.aotwoint, f2py=True)
        np.testing.assert_allclose(f_a, self.faref)
        np.testing.assert_allclose(f_b, self.fbref)

    def test_f_p(self):
        "Test total Fock, Python version"""
        dtot = self.daref + self.dbref
        ftot = fock(dtot, filename=self.aotwoint, f2py=False)

        fref = 0.5*(self.faref+self.fbref)
        np.testing.assert_allclose(ftot, fref)

    def test_f_f(self):
        "Test total Fock, Fortran version"""
        dtot = self.daref + self.dbref
        ftot = fock(dtot, filename=self.aotwoint, f2py=True)

        fref = 0.5*(self.faref + self.fbref)
        np.testing.assert_allclose(ftot, fref)

    def test_fs_p(self):
        "Test spin Fock, Python version"""
        dspin = self.daref - self.dbref
        fspin = fock(dspin, hfc=0, filename=self.aotwoint,  f2py=False)

        fref = 0.5*(self.faref - self.fbref)

        np.testing.assert_allclose(fspin, fref, atol=1e-8)

    def test_fs_f(self):
        "Test spin Fock, Python version"""
        dspin = self.daref - self.dbref
        fspin = fock(dspin, hfc=0, filename=self.aotwoint, f2py=True)

        fref = 0.5*(self.faref - self.fbref)

        np.testing.assert_allclose(fspin, fref, atol=1e-8)
