#!/usr/bin/env python
"""Module for getting stuff from Dalton two-electron integral file AOTWOINT"""

import numpy as np
from fortran_binary import FortranBinary as FB
from util.full import Matrix


def info(filename="AOTWOINT"):
    """Extract info block BASINFO on AOTWOINT"""
    _aotwoint = FB(filename)
    _aotwoint.find("BASINFO")
    rec = next(_aotwoint)
    fileinfo = rec.read(12,'i')
    retinfo = {
        "nsym":fileinfo[0],
        "nbas":fileinfo[1:9],
        "lbuf":fileinfo[9],
        "nibuf":fileinfo[10],
        "nbits":fileinfo[11]
        }
    return retinfo


def list_buffers(filename="AOTWOINT", label="BASTWOEL"):
    """ Return integral buffers in AOTWOINT"""
    _aofile = FB(filename)
    _aofile.find(label)

    for rec in _aofile:
        lbuf = (len(rec)-4) // 12

        buf = np.array(rec.read(lbuf,'d'))
        ibuf = np.array(rec.read(4*lbuf,'B')).reshape(lbuf, 4)
        length = rec.read(1,'i')[0]

        if length < 0: return
        yield buf[:length], ibuf[:length]

def list_integrals(*args, **kwargs):
    """ List two-electron spin-orbit integrals in file """

    for buf, ibuf in list_buffers(*args, **kwargs):
        for g, ig in zip(buf, ibuf):
            yield ig, g


def fockab(*Dab, **kwargs):
    """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

    assert len(Dab)>0

    f2py = kwargs.get('f2py', True)

    if f2py is True:
        try:
            import sirfck
        except(ImportError):#pragma: no cover
            f2py = False
            print("Warning: non-existent sirfck.so wanted - reverting to python")

    if f2py:
        Fab = fock_builder_f(Dab, **kwargs)
    else:
        Fab = fock_builder_py(Dab, **kwargs)
    return Fab

def fock_builder_f(Dab, **kwargs):
    import sirfck

    filename = kwargs.get('filename', 'AOTWOINT')
    hfc = kwargs.get('hfc', 1)
    hfx = kwargs.get('hfx', 1)

    Dshape = Dab[0][0].shape
    fdim = (*Dshape, len(Dab))
    Js = Matrix(fdim)
    Kas = Matrix(fdim)
    Kbs = Matrix(fdim)
    Das = Matrix(fdim)
    Dbs = Matrix(fdim)
    for i, (Da, Db) in enumerate(Dab):
        Das[:, :, i] = Da
        Dbs[:, :, i] = Db
    for buf, ibuf in list_buffers(filename):
        Js, Kas, Kbs = sirfck.fckab(
            Js, Kas, Kbs, Das, Dbs, buf, ibuf.T
            )

    Fas = (hfc*Js[:, :, i] - hfx*Kas[:, :, i] for i in range(len(Dab)))
    Fbs = (hfc*Js[:, :, i] - hfx*Kbs[:, :, i] for i in range(len(Dab)))

    Fab = tuple( (Fa, Fb) for Fa, Fb in zip(Fas, Fbs))
    return Fab

def fock_builder_py(Dab, **kwargs):

    filename = kwargs.get('filename', 'AOTWOINT')
    hfc = kwargs.get('hfc', 1)
    hfx = kwargs.get('hfx', 1)

    Ds = [Da + Db for Da, Db in Dab]
    Js = [Matrix(D.shape) for D in Ds]
    Ks = [(Matrix(Da.shape), Matrix(Db.shape)) for Da, Db in Dab]
    for ig, g in list_integrals(filename):
        p, q, r, s = ig
        s, r, q, p = (p-1, q-1, r-1, s-1)
        if ( p == q ): g *= .5
        if ( r == s ): g *= .5
        if ( p == r and q == s ): g *= .5

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

def fock(D, filename="AOTWOINT", hfc=1, hfx=1, f2py=True):
    """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

    if f2py is True:
        try:
            import sirfck
        except(ImportError):#pragma: no cover
            f2py = False
            print("Warning: non-existent sirfck.so wanted - reverting to python")

    J = Matrix(D.shape)
    K = Matrix(D.shape)

    if f2py:
        for buf, ibuf in list_buffers(filename):
            J, K = sirfck.fck(
                J, K,  D, D, buf, ibuf.T
                )
    else:
        for ig, g in list_integrals(filename):
            p, q, r, s = ig
            s, r, q, p = (p-1, q-1, r-1, s-1)
            if ( p == q ): g *= .5
            if ( r == s ): g *= .5
            if ( p == r and q == s ): g *= .5

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



def main():

    import argparse 

    parser = argparse.ArgumentParser()
    parser.add_argument('--file', default='/tmp/AOTWOINT', help='Two-electron integral file')
    parser.add_argument( '--list', action='store_true', help='List integrals on file')

    args = parser.parse_args()

    if args.list:
        print("List integrals")
        for ig, g in list_integrals(args.file):
            print(ig, g)

if __name__ == "__main__":#pragma: no cover
    import sys
    sys.exit(main())
