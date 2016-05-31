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

        norb = [9, 4, 4, 2, 4, 2, 2, 0]
        self.iorb = [sum(norb[:i]) for i in range(8)]
        self.occa = [3, 2, 2, 0, 2, 0, 0, 0]
        self.occb = [3, 1, 2, 0, 2, 0, 0, 0]

        nbast = sum(norb)
        self.nbast = nbast
        self.da = self.get_dens(self.occa)
        self.db = self.get_dens(self.occb)

    def get_dens(self, ndim):
        dens = matrix((self.nbast, self.nbast))
        for ni, offset in zip(ndim, self.iorb):
            if ni > 0:
                cmoi = self.cmo[:, offset: offset + ni]
                dens += cmoi*cmoi.T

        return dens


class TestSpinOrbitCla(TestSpinOrbitCl):


    def setUp(self):
        TestSpinOrbitCl.setUp(self)

        nish = [3, 1, 1, 0, 1, 0, 0, 0]
        nash = [0, 1, 1, 0, 1, 0, 0, 0]
        self.di = 2*self.get_dens(nish)
        self.dv = 2*self.get_dens(nash)


    def test_fc(self):
        fc = twoso.fock(self.di, 'z' , filename=self.ao2soint)
        fcmo = self.cmo.T*fc*self.cmo
        np.testing.assert_almost_equal(fcmo[10, 14], -20.80725730)

    @unittest.skip('holdon')
    def test_fab(self):
        fa, fb = twoso.fockab(self.da, self.db, 'z', filename=self.ao2soint)
        fcmo = self.cmo.T*(fa-fb)*self.cmo
        np.testing.assert_almost_equal(fcmo[10, 14], -20.80725730)


class TestSpinOrbitClb(TestSpinOrbitCl):

    def setUp(self):
        TestSpinOrbitCl.setUp(self)

        nish = [3, 1, 1, 0, 2, 0, 0, 0]
        nash = [0, 1, 1, 0, 0, 0, 0, 0]
        self.di = 2*self.get_dens(nish)

    def test_fc(self):
        fc = twoso.fock(self.di, 'z' , filename=self.ao2soint)
        fcmo = self.cmo.T*fc*self.cmo
        np.testing.assert_almost_equal(fcmo[10, 14], -21.29974013)


if __name__ == "__main__":
    unittest.main()
