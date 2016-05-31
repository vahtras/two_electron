import unittest
import mock
import sys
import os
import numpy as np
from util.full import init, matrix
from . import two
from two import twoso

class TestSpinOrbitCl(unittest.TestCase):
    def setUp(self):
        n, _ = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.ao2soint = os.path.join(suppdir, "AO2SOINT")
        self.cmo = np.loadtxt(os.path.join(suppdir, "cmo")).view(matrix)

        norb = np.array([9, 4, 4, 2, 4, 2, 2, 0])
        self.iorb = np.array([sum(norb[:i]) for i in range(8)])
        self.occa = np.array([3, 2, 2, 0, 2, 0, 0, 0])
        self.occb = np.array([3, 1, 2, 0, 2, 0, 0, 0])

        nbast = sum(norb)
        self.nbast = nbast
        self.da = self.get_dens(self.occa, self.iorb)
        self.db = self.get_dens(self.occb, self.iorb)

    def get_dens(self, ndim, idim, nocc=None):
        if nocc is None:
            nocc_ = [1]*len(ndim)
        else:
            nocc_ = nocc

        dens = matrix((self.nbast, self.nbast))
        for ni, no, offset in zip(ndim, nocc_, idim):
            if ni > 0:
                cmoi = self.cmo[:, offset: offset + ni]
                dens += no*cmoi*cmoi.T

        return dens


class TestSpinOrbitCla(TestSpinOrbitCl):


    def setUp(self):
        TestSpinOrbitCl.setUp(self)

        nish = np.array([3, 1, 1, 0, 1, 0, 0, 0])
        nash = np.array([0, 1, 1, 0, 1, 0, 0, 0])
        self.di = 2*self.get_dens(nish, self.iorb)
        self.dv = self.get_dens(
            nash,
            self.iorb + nish,
            nocc=[0, 1, 2, 0, 2, 0, 0, 0]
            )

    def test_internal(self):
        nish = [3, 1, 1, 0, 1, 0, 0, 0]
        nocc = [2, 2, 2, 0, 2, 0, 0, 0]
        np.testing.assert_allclose(
            2*self.get_dens(nish, self.iorb), 
            self.get_dens(nish, self.iorb, nocc)
            )

    def test_internal_densities(self):
        np.testing.assert_allclose(self.da+self.db, self.di+self.dv)


    def test_fc(self):
        fc = twoso.fock(self.di, 'z' , filename=self.ao2soint)
        fcmo = self.cmo.T*fc*self.cmo
        np.testing.assert_almost_equal(fcmo[10, 14], -20.80725730)

    def test_fab(self):
        fa, fb = twoso.fockab(self.da, self.db, 'z', filename=self.ao2soint)
        fc = twoso.fock(self.di, 'z', filename=self.ao2soint)
        fv = twoso.fock(self.dv, 'z', filename=self.ao2soint)
        np.testing.assert_almost_equal(0.5*(fa - fb), fc + fv)


class TestSpinOrbitClb(TestSpinOrbitCl):

    def setUp(self):
        TestSpinOrbitCl.setUp(self)

        nish = [3, 1, 1, 0, 2, 0, 0, 0]
        nash = [0, 1, 1, 0, 0, 0, 0, 0]
        self.di = 2*self.get_dens(nish, self.iorb)

    def test_fc(self):
        fc = twoso.fock(self.di, 'z' , filename=self.ao2soint)
        fcmo = self.cmo.T*fc*self.cmo
        np.testing.assert_almost_equal(fcmo[10, 14], -21.29974013)


if __name__ == "__main__":
    unittest.main()
