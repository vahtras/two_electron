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
        assert ig  == (1, 1, 1, 1)


class TestH2O:
    def setup_method(self):
        suppdir  = pathlib.Path(__file__).with_suffix(".d") / "H2O"

        self.aotwoint = suppdir / "hf_H2O_ccpVDZ.AOTWOINT"
        self.db = suppdir / "aotwoint.db"
        self.reader = two.eri.Reader(self.aotwoint)
        self.freader = two.eri.FReader(self.aotwoint)
        self.sqlreader = two.eri.SQLReader(self.aotwoint, db=self.db)

        self.d = np.loadtxt(suppdir / 'dcao').reshape((24, 24))
        self.f = np.loadtxt(suppdir / 'fcao').reshape((24, 24))

    def test_number_of_integrals_orig(self):
        assert len(list(self.reader.list_integrals())) ==  11412

    def test_number_of_integrals_sql(self):
        # self.sqlreader.insert_integrals()
        assert len(list(self.sqlreader.list_integrals())) ==  11412

    def test_dens_fock_py(self):
        fock = self.reader.fock(self.d)
        np.testing.assert_almost_equal(fock, self.f)

    def test_dens_fock_f(self):
        fock = self.freader.fock(self.d)
        np.testing.assert_almost_equal(fock, self.f)

    @mark.skip
    def test_dens_fock_sql(self):
        fock = self.sqlreader.fock(self.d)
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
