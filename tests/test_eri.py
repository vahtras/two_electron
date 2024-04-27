import pathlib
import sys
from unittest.mock import patch

from pytest import approx, mark
import numpy as np

import two.eri

class TestERI:
    def setup_method(self):
        suppdir  = pathlib.Path(__file__).with_suffix(".d")
        self.aotwoint = suppdir / "AOTWOINT"
        self.aotwoint_db = suppdir / "aotwoint.db"
        self.reader = two.eri.Reader(self.aotwoint)
        self.freader = two.eri.FReader(self.aotwoint)
        self.sqlreader = two.eri.SQLReader(self.aotwoint, db=None)
        self.sqlreader.insert_integrals()

    def test_basinfo(self):
        info = self.reader.basinfo()
        assert info["nsym"] == 1
        assert info["nbas"] == (7, 0, 0, 0, 0, 0, 0, 0)
        assert info["lbuf"] == 600
        assert info["nibuf"] == 1
        assert info["nbits"] == 8

    @mark.parametrize(
        'reader',
        ["reader", "freader", "sqlreader"]
    ) 
    def test_first_integral(self, reader):
        for ig, g in getattr(self, reader).list_integrals():
            break
        assert g == approx(4.78506540471)
        if reader == 'sqlreader':
            assert ig  == (0, 0, 0, 0)
        else:
            assert ig  == (1, 1, 1, 1)

    def test_coulomb(self):
        dens = np.eye(7)
        f1 = self.reader.fock(dens, hfx=0)
        f2 = self.sqlreader.fock(dens, hfx=0)
        np.testing.assert_almost_equal(f1, f2)

    def test_exchange(self):
        dens = np.eye(7)
        f1 = self.reader.fock(dens, hfc=0, hfx=-2)
        f2 = self.sqlreader.fock(dens, hfc=0, hfx=-2)
        np.testing.assert_almost_equal(f1, f2)

    def test_total_fock(self):
        dens = np.eye(7)
        f1 = self.reader.fock(dens)
        f2 = self.sqlreader.fock(dens)
        np.testing.assert_almost_equal(f1, f2)



class TestH2O:
    def setup_method(self):
        suppdir  = pathlib.Path(__file__).with_suffix(".d") / "H2O"

        self.aotwoint = suppdir / "hf_H2O_ccpVDZ.AOTWOINT"
        self.db = suppdir / "aotwoint.db"
        self.reader = two.eri.Reader(self.aotwoint)
        self.freader = two.eri.FReader(self.aotwoint)
        self.sqlreader = two.eri.SQLReader(self.aotwoint, db=self.db)
        self.sqlreader.insert_integrals()

        self.d = np.loadtxt(suppdir / 'dcao').reshape((24, 24))
        self.f = np.loadtxt(suppdir / 'fcao').reshape((24, 24))

    @mark.parametrize(
        'reader',
        ["reader", "sqlreader"]
    ) 
    def test_number_of_integrals(self, reader):
        assert len(list(getattr(self,reader).list_integrals())) ==  11412


    @mark.parametrize(
        'reader',
        ["reader", "freader"]#, "sqlreader"]
    ) 
    def test_dens_fock(self, reader):
        fock = getattr(self, reader).fock(self.d)
        np.testing.assert_almost_equal(fock, self.f)

    @patch('two.eri.Reader.list_integrals')
    def test_arg_int(self, mock_list):
        sys.argv[1:] = ['/dev/null/AOTWOINT', '--list']
        mock_list.return_value=[((0,0,0,0), 3.14)]
        two.eri.main()
        mock_list.assert_called_once_with()

    # @mark.skip
    def test_insert_integrals(self):
        new = two.eri.SQLReader(self.aotwoint, db='/tmp/foo.db')
        new.insert_integrals()
        assert len(list(new.list_integrals())) ==  11412

    def test_insert_density(self):
        new = two.eri.SQLReader(self.aotwoint, db='/tmp/foo.db')
        new.insert_density(self.d)
        assert len(list(new.list_density())) ==  24*24
