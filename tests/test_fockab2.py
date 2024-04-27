import pathlib

import numpy as np
import two.eri

from pytest import mark


@mark.parametrize(
    'reader',
    [
        ('reader'),
        ('freader'),
    ]
)
class TestTwo:
    """Two-fock test suite"""

    def setup_method(self):
        """Setup supporting directory"""
        self.suppdir = pathlib.Path(__file__).with_suffix(".d")
        self.aotwoint = self.suppdir/"AOTWOINT"
        self.reader = two.eri.Reader(self.aotwoint)
        self.freader = two.eri.FReader(self.aotwoint)

        self.faref = np.loadtxt(self.suppdir/'fa.ref')
        self.fbref = np.loadtxt(self.suppdir/'fb.ref')

        self.daref = np.loadtxt(self.suppdir/'da.ref')
        self.dbref = np.loadtxt(self.suppdir/'db.ref')

    def test_fab(self, reader):
        """Test alpha and beta Fock matrix"""
        d_a, d_b = self.daref, self.dbref
        (f_a, f_b), = getattr(self, reader).fockab((d_a, d_b))
        np.testing.assert_allclose(f_a, self.faref)
        np.testing.assert_allclose(f_b, self.fbref)

    def test_f(self, reader):
        "Test total Fock, Python version"""
        dtot = self.daref + self.dbref
        ftot = getattr(self, reader).fock(dtot)

        fref = 0.5*(self.faref+self.fbref)
        np.testing.assert_allclose(ftot, fref)


    def test_fs(self, reader):
        "Test spin Fock, Python version"""
        dspin = self.daref - self.dbref
        fspin = getattr(self, reader).fock(dspin, hfc=0)

        fref = 0.5*(self.faref - self.fbref)

        np.testing.assert_allclose(fspin, fref, atol=1e-8)
