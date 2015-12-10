import  unittest
import os
import subprocess
import numpy
from ..util.full import matrix
from .. import two
import mock

class TwoTest(unittest.TestCase):

    def setUp(self):
        n, e = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.aotwoint = os.path.join(suppdir, "AOTWOINT")

    def test_basinfo(self):
        info = two.info(self.aotwoint)
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

class TwoIntTest(unittest.TestCase):

    def setUp(self):
        n, e = os.path.splitext(__file__)
        suppdir = n + ".d"
        self.twoint = two.TwoInt(os.path.join(suppdir, "AOTWOINT"))
        self.info = self.twoint.info()

    def test_basinfo_nsym(self):
        self.assertEqual(self.info["nsym"], 1)

    def test_basinfo_nbas(self):
        numpy.testing.assert_equal(self.info["nbas"], [7,0,0,0,0,0,0,0])

    def test_basinfo_lbuf(self):
        self.assertEqual(self.info["lbuf"], 600)

    def test_basinfo_nibuf(self):
        self.assertEqual(self.info["nibuf"], 1)

    def test_basinfo_nbits(self):
        self.assertEqual(self.info["nbits"], 8)
    
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
            #subprocess.call(args)

    @unittest.skip('generate integral files to run')
    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 11412)

    def setUp(self):
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((24, 24))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((24, 24))

    @unittest.skip('generate integral files to run')
    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=False), 
            self.f
            )

    @unittest.skip('generate integral files to run')
    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=True), 
            self.f
        )

class TestTwoIntH2o(unittest.TestCase):

    def setUp(self):
        return
        root, ext = n, e = os.path.splitext(__file__)
        self.tmpdir = os.path.join(root + ".d", "H2O")
        self.aotwoint = two.TwoInt(os.path.join(self.tmpdir, "hf_H2O_ccpVDZ.AOTWOINT"))
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((24, 24))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((24, 24))

    @unittest.skip('generate integral files to run')
    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(
            self.aotwoint.fock(self.d, f2py=False), 
            self.f
            )

    @unittest.skip('generate integral files to run')
    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            self.aotwoint.fock(self.d, f2py=True), 
            self.f
        )


if __name__ == "__main__":
    unittest.main()

    
