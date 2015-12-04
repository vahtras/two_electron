import  unittest
import os
import subprocess
import numpy
from util.full import matrix
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
    
    @mock.patch('two_electron.two.TwoInt.list_buffers')
    def test_list(self, mock_buffers):
        mock_buffers.return_value = ([], [])
        d = numpy.ones((7, 7))
        f = self.twoint.fock(d, hfx=0)
        numpy.testing.assert_allclose(f, d)


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
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=False), 
            self.f
            )

    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=True), 
            self.f
        )

class TestTwoIntH2o(unittest.TestCase):

    def setUp(self):
        root, ext = n, e = os.path.splitext(__file__)
        self.tmpdir = os.path.join(root + ".d", "H2O")
        self.aotwoint = two.TwoInt(os.path.join(self.tmpdir, "hf_H2O_ccpVDZ.AOTWOINT"))
        self.d = numpy.loadtxt(os.path.join(self.tmpdir, 'dcao')).view(matrix).reshape((24, 24))
        self.f = numpy.loadtxt(os.path.join(self.tmpdir, 'fcao')).view(matrix).reshape((24, 24))

    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(
            self.aotwoint.fock(self.d, f2py=False), 
            self.f
            )

    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            self.aotwoint.fock(self.d, f2py=True), 
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

    @unittest.skip('long test')
    def test_number_of_integrals(self):
        self.assertEqual(len(list(two.list_integrals(self.aotwoint))), 972549)

    @unittest.skip('long test')
    def test_dens_fock(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=False), 
            self.f
            )

    @unittest.skip('long test')
    def test_dens_fock_f2py(self):
        numpy.testing.assert_almost_equal(
            two.fock(self.d, filename=self.aotwoint, f2py=True), 
            self.f
        )

class TestTransform(unittest.TestCase):

    @mock.patch('two_electron.two.matrix')
    @mock.patch('two_electron.two.list_integrals')
    def test_matrix_call(self, mock_list_integrals, mock_matrix):
	d1=mock.Mock(spec=['shape'])
	d1.shape=(30,60)
	mmaa = two.semitransform(d1, d1)
	mock_matrix.assert_called_with((30, 30, 60, 60))
	

    @mock.patch('two_electron.two.list_integrals')
    def test_list_integrals_call(self, mock_list_integrals):
	d1=mock.Mock(spec=['shape'])
	d1.shape=(30,60)
	two.semitransform(d1, d1, file='integral_file')
	mock_list_integrals.assert_called_once_with('integral_file')

    @mock.patch('two_electron.two.list_integrals')
    def test_list_integrals_one(self, mock_list_integrals):
	d1 = matrix((1, 1))
	d1[0, 0] = 2
	mock_list_integrals.return_value = [((1, 1, 1, 1), 3.0)]
	mmaa = two.semitransform(d1, d1, file='integral_file')
	numpy.testing.assert_allclose([[[[12]]]], mmaa)

    @mock.patch('two_electron.two.list_integrals')
    def test_list_integrals_two(self, mock_list_integrals):
	d1 = matrix((1, 2))
	d1[0, 0] = d1[0, 1] = 1
	mock_list_integrals.return_value = [
	    ((1, 1, 1, 1), 0.0),
	    ((2, 1, 1, 1), 1.0),
	    ((2, 1, 2, 1), 2.0),
	    ((2, 2, 1, 1), 2.0),
	    ((2, 2, 2, 1), 3.0),
	    ((2, 2, 2, 2), 4.0)
	    ]
	mmaa = two.semitransform(d1, d1, file='integral_file')
	numpy.testing.assert_allclose([[[[4, 8], [8, 12]]]], mmaa)

if __name__ == "__main__":
    unittest.main()

    
