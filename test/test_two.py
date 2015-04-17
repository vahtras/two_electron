import  unittest
import os
import subprocess
import numpy
from util.full import matrix
from .. import two

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
            #os.system('dalton -get AOTWOINT hf H2O_ccpVDZ')
            args = ['dalton', '-get', 'AOTWOINT', cls.dal, cls.mol]
            subprocess.call(args)

    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 11412)

    def setUp(self):
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((24, 24))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((24, 24))

    def test_dens_fock(self):
	print self.d
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=False), 
            self.f
            )

    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=True), 
            self.f
        )

class TestAcetaldehyde(TestBase):


    @classmethod
    def setUpClass(cls):
        super(TestAcetaldehyde, cls).setUpClass()
        cls.subdir = "acetaldehyde"
        cls.dal = "hf"
        cls.mol = "acetaldehyde"
        cls.filename = "%s_%s.AOTWOINT" % (cls.dal, cls.mol)
        cls.tmpdir = os.path.join(cls.base_dir, cls.subdir)
        cls.aotwoint = os.path.join(cls.tmpdir, cls.filename)
        if not os.path.exists(cls.aotwoint):
            os.chdir(cls.tmpdir)
            #os.system('dalton -get AOTWOINT %s %s' % (cls.dal, cls.mol))
            args = ['dalton', '-get', 'AOTWOINT', cls.dal, cls.mol]
            subprocess.call(args)

    def setUp(self):
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((146, 146))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((146, 146))

    @unittest.skip('long test')
    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 28346779)

    @unittest.skip('long test')
    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=False), 
            self.f
            )

    @unittest.skip('segfaults')
    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=True), 
            self.f
        )


class TestAcetaldehydeSmall(TestBase):

    @classmethod
    def setUpClass(cls):
        super(TestAcetaldehydeSmall, cls).setUpClass()
        cls.subdir = "acetaldehyde_small"
        cls.dal = "hf"
        cls.mol = "acetaldehyde"
        cls.filename = "%s_%s.AOTWOINT" % (cls.dal, cls.mol)
        cls.tmpdir = os.path.join(cls.base_dir, cls.subdir)
        cls.aotwoint = os.path.join(cls.tmpdir, cls.filename)
        if not os.path.exists(cls.aotwoint):
            os.chdir(cls.tmpdir)
            #os.system('dalton -get AOTWOINT %s %s' % (cls.dal, cls.mol))
            args = ['dalton', '-get', 'AOTWOINT', cls.dal, cls.mol]
            subprocess.call(args)

    def setUp(self):
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((62, 62))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((62, 62))

    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 972549)

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

if __name__ == "__main__":
    unittest.main()

    
