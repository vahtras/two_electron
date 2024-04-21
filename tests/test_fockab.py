"""Some docstring"""
import os
import pathlib

import numpy as np
import two.eri
import pytest


class TestTwo:
    """Two-fock test suite"""

    def setup_method(self):
        """Setup supporting directory"""
        name, _ = os.path.splitext(__file__)
        self.suppdir = pathlib.Path(name + ".d")
        self.aotwoint = self.suppdir/"AOTWOINT"
        self.reader = two.eri.Reader(self.aotwoint)
        self.freader = two.eri.FReader(self.aotwoint)

        self.daref = np.array([
            [1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            ])

        self.dbref = np.array([
            [1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
            ])

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

    @pytest.mark.parametrize("reader", ["reader", "freader"])
    def test_fab(self, reader):
        """Test alpha and beta Fock matrix"""
        (f_a, f_b), = getattr(self, reader).fockab((self.daref, self.dbref))
        np.testing.assert_allclose(f_a, self.faref)
        np.testing.assert_allclose(f_b, self.fbref)

    @pytest.mark.parametrize("reader", ["reader", "freader"])
    def test_f(self, reader):
        "Test total Fock, Python version"""
        dtot = self.daref + self.dbref
        ftot = getattr(self, reader).fock(dtot)

        fref = 0.5*(self.faref+self.fbref)
        np.testing.assert_allclose(ftot, fref)

    @pytest.mark.parametrize("reader", ["reader", "freader"])
    def test_fs(self, reader):
        "Test spin Fock, Python version"""
        dspin = self.daref - self.dbref
        fspin = getattr(self, reader).fock(dspin, hfc=0)

        fref = 0.5*(self.faref-self.fbref)
        np.testing.assert_allclose(fspin, fref, atol=1e-8)
