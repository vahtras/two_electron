import unittest
try:
    import mock
except ImportError:
    import unittest.mock as mock
import sys
import os
import numpy as np
from util.full import init
from . import two
from two import twoso

class TestSpinOrbit(unittest.TestCase):


    def setUp(self):
        n, _ = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.ao2soint = os.path.join(suppdir, "AO2SOINT")

        self.dc = init([
           [ 2.10662148,    -0.44832042,     0.00000000,     0.10863721,     0.00000000,     -0.02702924,    -0.02702924],
           [-0.44832042,     1.97509195,    -0.00000000,    -0.61381008,     0.00000000,     -0.03765992,    -0.03765992],
           [ 0.00000000,    -0.00000000,     0.73437339,     0.00000000,     0.00000000,      0.54049192,    -0.54049192],
           [ 0.10863721,    -0.61381008,     0.00000000,     1.22749198,     0.00000000,      0.47809924,     0.47809924],
           [ 0.00000000,     0.00000000,     0.00000000,     0.00000000,     2.00000000,      0.00000000,     0.00000000],
           [-0.02702924,    -0.03765992,     0.54049192,     0.47809924,     0.00000000,      0.60856871,    -0.18702525],
           [-0.02702924,    -0.03765992,    -0.54049192,     0.47809924,     0.00000000,     -0.18702525,     0.60856871] 
        ])


        self.fcref = -init([
            [ 0.00000000,     0.00000000,     0.37936213,    -0.00000000,     0.00000000,      0.00253951,    -0.00253951],
            [-0.00000000,     0.00000000,     0.10878631,    -0.00000000,     0.00000000,     -0.01699688,     0.01699688],
            [-0.37936213,    -0.10878631,     0.00000000,   -11.60434639,     0.00000000,     -0.98565692,    -0.98565692],
            [ 0.00000000,     0.00000000,    11.60434639,     0.00000000,     0.00000000,      1.20385366,    -1.20385366],
            [ 0.00000000,     0.00000000,     0.00000000,     0.00000000,     0.00000000,      0.00000000,     0.00000000],
            [-0.00253951,     0.01699688,     0.98565692,    -1.20385366,     0.00000000,      0.00000000,    -0.36840521],
            [ 0.00253951,    -0.01699688,     0.98565692,     1.20385366,     0.00000000,      0.36840521,     0.00000000] 
        ])

    def test_first_integral(self):
        for c, ig, g in twoso.list_integrals(self.ao2soint):
            break
        self.assertEqual(c, 0)
        self.assertAlmostEqual(g, 2.88719908251)

    def test_fc(self):
        fc = twoso.fock(self.dc, 'z' , filename=self.ao2soint)
        np.testing.assert_almost_equal(fc, self.fcref)

    def test_fab(self):
        da = db = 0.5*self.dc
        fa, fb  = twoso.fockab(da, db, 'z' , filename=self.ao2soint)
        fc = 0.5*(fa - fb)
        np.testing.assert_almost_equal(fc, self.fcref)

    @mock.patch('two.twoso.list_integrals')
    def test_arg_int(self, mock_list):
        print(sys.argv)
        sys.argv[1:] = ['/dev/null/AO2SOINT', '-l']
        mock_list.return_value=[('x', (0,0,0,0), 3.14)]
        twoso.main()
        mock_list.assert_called_once_with('/dev/null/AO2SOINT')

if __name__ == "__main__":
    unittest.main()
