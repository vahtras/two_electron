import os
import subprocess
import sys
import unittest

import numpy
from util.full import matrix

import two.core

try:
    import mock
except ImportError:
    import unittest.mock as mock

class TwoTest(unittest.TestCase):

    def setUp(self):
        n, e = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.aotwoint = os.path.join(suppdir, "AOTWOINT")

    def test_basinfo(self):
        info = two.core.info(self.aotwoint)
        self.assertEqual(info["nsym"], 1)
        numpy.testing.assert_equal(info["nbas"], [7,0,0,0,0,0,0,0])
        self.assertEqual(info["lbuf"], 600)
        self.assertEqual(info["nibuf"], 1)
        self.assertEqual(info["nbits"], 8)


    def test_first_integral(self):
        for ig, g in two.list_integrals(self.aotwoint):
            break
        self.assertAlmostEqual(g, 4.78506540471)
        numpy.testing.assert_equal(ig, [1,1,1,1])

    
class TestBase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        root, ext = n, e = os.path.splitext(__file__)
        cls.base_dir = root + ".d"
        
class TestH2O(TestBase):

    @classmethod
    def setUpClass(cls):
        super(TestH2O, cls).setUpClass()
        cls.subdir = "H2O"
        cls.dal = "hf"
        cls.mol = "H2O_ccpVDZ"
        cls.filename = "%s_%s.AOTWOINT" % (cls.dal, cls.mol)
        cls.tmpdir = os.path.join(cls.base_dir, cls.subdir)
        cls.aotwoint = os.path.join(cls.tmpdir, cls.filename)
        if not os.path.exists(cls.aotwoint):
            os.chdir(cls.tmpdir)
            args = ['dalton', '-get', 'AOTWOINT', cls.dal, cls.mol]
            subprocess.call(args)

    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 11412)

    def setUp(self):
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((24, 24))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((24, 24))

    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=False), 
            self.f
            )

    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=True), 
            self.f
        )

    @mock.patch('two.core.list_integrals')
    def test_arg_int(self, mock_list):
        print(sys.argv)
        sys.argv[1:] = ['/dev/null/AOTWOINT', '--list']
        mock_list.return_value=[((0,0,0,0), 3.14)]
        two.main()
        mock_list.assert_called_once_with('/dev/null/AOTWOINT')


if __name__ == "__main__":
    unittest.main()

    
