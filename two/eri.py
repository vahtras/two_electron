import sys

from fortran_binary import FortranBinary
import numpy as np
try:
    import sirfck
except ImportError:
    pass

class Reader:
    def __init__(self, filename):
        self.filename = filename

    def basinfo(self):
        """
        Extract info block BASINFO from AOTWOINT file.
        """
        with FortranBinary(self.filename) as f:
            f.find("BASINFO")
            rec = next(f)

        fileinfo = rec.read(12, 'i')
        retinfo = {
            "nsym": fileinfo[0],
            "nbas": fileinfo[1:9],
            "lbuf": fileinfo[9],
            "nibuf": fileinfo[10],
            "nbits": fileinfo[11],
        }
        return retinfo

    def list_integrals(self):
        """
        List two-electron spin-orbit integrals in file
        """
        for buf, ibuf in self.list_buffers():
            for g, ig in zip(buf, ibuf):
                yield ig, g

    def list_buffers(self, label="BASTWOEL"):
        """
        Return integral buffers in AOTWOINT
        """
        _aofile = FortranBinary(self.filename)
        _aofile.find(label)

        for rec in _aofile:
            lbuf = (len(rec)-4) // 12

            buf = np.array(rec.read(lbuf, 'd'))
            ibuf = np.array(rec.read(4*lbuf, 'B')).reshape(lbuf, 4)
            length = rec.read(1, 'i')[0]

            if length < 0:
                return
            yield buf[:length], ibuf[:length]

    def fock(self, D, hfc=1, hfx=1, f2py=True):
        """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

        if 'sirfck' in sys.modules:
            pass
        else:
            import warnings
            warnings.warn(
                "Warning: non-existent sirfck.so wanted"
                "- reverting to python"
            )
            f2py = False

        J = np.zeros(D.shape)
        K = np.zeros(D.shape)

        if f2py:
            try:
                import sirfck
            except ImportError:
                print("Warning: sirfck.so not found - reverting to python")
                f2py = False

            for buf, ibuf in self.list_buffers(self.filename):
                J, K = sirfck.fck(
                    J, K,  D, D, buf, ibuf.T
                    )
        else:
            for ig, g in self.list_integrals():
                p, q, r, s = ig
                s, r, q, p = (p-1, q-1, r-1, s-1)
                if (p == q):
                    g *= .5
                if (r == s):
                    g *= .5
                if (p == r and q == s):
                    g *= .5

                Jadd = g * (D[r, s] + D[s, r])
                J[p, q] += Jadd
                J[q, p] += Jadd
                Jadd = g * (D[p, q] + D[q, p])
                J[r, s] += Jadd
                J[s, r] += Jadd
                K[p, s] += g*D[r, q]
                K[p, r] += g*D[s, q]
                K[q, s] += g*D[r, p]
                K[q, r] += g*D[s, p]
                K[r, q] += g*D[p, s]
                K[s, q] += g*D[p, r]
                K[r, p] += g*D[q, s]
                K[s, p] += g*D[q, r]

        return hfc*J - 0.5*hfx*K

    def fockab(self, *Dab, **kwargs):
        """
        Routine for building alpha and beta fock matrix by reading AOTWOINT
        """

        assert len(Dab) > 0

        f2py = kwargs.get('f2py', True)

        if f2py:
            if 'sirfck' in sys.modules:
                pass
            else:
                print(
                    "Warning: non-existent sirfck.so wanted"
                    "- reverting to python"
                )
                f2py = False

        if f2py:
            Fab = self.fock_builder_f(Dab, **kwargs)
        else:
            Fab = self.fock_builder_py(Dab, **kwargs)
        return Fab

    def fock_builder_f(self, Dab, **kwargs):

        filename = kwargs.get('filename', 'AOTWOINT')
        hfc = kwargs.get('hfc', 1)
        hfx = kwargs.get('hfx', 1)

        Dshape = Dab[0][0].shape
        fdim = (*Dshape, len(Dab))
        Js = np.zeros(fdim)
        Kas = np.zeros(fdim)
        Kbs = np.zeros(fdim)
        Das = np.zeros(fdim)
        Dbs = np.zeros(fdim)
        for i, (Da, Db) in enumerate(Dab):
            Das[:, :, i] = Da
            Dbs[:, :, i] = Db
        for buf, ibuf in self.list_buffers(filename):
            Js, Kas, Kbs = sirfck.fckab(
                Js, Kas, Kbs, Das, Dbs, buf, ibuf.T
                )

        Fas = (hfc*Js[:, :, i] - hfx*Kas[:, :, i] for i in range(len(Dab)))
        Fbs = (hfc*Js[:, :, i] - hfx*Kbs[:, :, i] for i in range(len(Dab)))

        Fab = tuple((Fa, Fb) for Fa, Fb in zip(Fas, Fbs))
        return Fab


    def fock_builder_py(self, Dab, **kwargs):

        hfc = kwargs.get('hfc', 1)
        hfx = kwargs.get('hfx', 1)

        Ds = [Da + Db for Da, Db in Dab]
        Js = [np.zeros(D.shape) for D in Ds]
        Ks = [(np.zeros(Da.shape), np.zeros(Db.shape)) for Da, Db in Dab]
        for ig, g in self.list_integrals():
            p, q, r, s = ig
            s, r, q, p = (p-1, q-1, r-1, s-1)
            if (p == q):
                g *= .5
            if (r == s):
                g *= .5
            if (p == r and q == s):
                g *= .5

            for D, (Da, Db), J, (Ka, Kb) in zip(Ds, Dab, Js, Ks):
                Jadd = g * (D[r, s] + D[s, r])
                J[p, q] += Jadd
                J[q, p] += Jadd
                Jadd = g * (D[p, q] + D[q, p])
                J[r, s] += Jadd
                J[s, r] += Jadd
                Ka[p, s] += g*Da[r, q]
                Ka[p, r] += g*Da[s, q]
                Ka[q, s] += g*Da[r, p]
                Ka[q, r] += g*Da[s, p]
                Ka[r, q] += g*Da[p, s]
                Ka[s, q] += g*Da[p, r]
                Ka[r, p] += g*Da[q, s]
                Ka[s, p] += g*Da[q, r]
                Kb[p, s] += g*Db[r, q]
                Kb[p, r] += g*Db[s, q]
                Kb[q, s] += g*Db[r, p]
                Kb[q, r] += g*Db[s, p]
                Kb[r, q] += g*Db[p, s]
                Kb[s, q] += g*Db[p, r]
                Kb[r, p] += g*Db[q, s]
                Kb[s, p] += g*Db[q, r]

        Fas = (hfc*J-hfx*Ka for J, (Ka, _) in zip(Js, Ks))
        Fbs = (hfc*J-hfx*Kb for J, (_, Kb) in zip(Js, Ks))

        Fab = tuple((Fa, Fb) for Fa, Fb in zip(Fas, Fbs))
        return Fab

def main():

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'file',
        help='Two-electron integral file'
    )
    parser.add_argument(
        '-l', '--list',
        action='store_true',
        help='List integrals on file'
    )

    args = parser.parse_args()

    if args.list:
        print("List integrals")
        reader = Reader(args.file)
        for ig, g in reader.list_integrals():
            print(ig, g)


if __name__ == "__main__":
    sys.exit(main())
