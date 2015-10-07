#!/usr/bin/env python
"""Module for getting stuff from Dalton two-electron integral file AOTWOINT"""

import numpy as np
from util.unformatted import FortranBinary as FB
from util.full import matrix

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
    _aofile = FB(filename, label=label)

    for rec in _aofile:
        lbuf = (_aofile.reclen-4)/12

        buf = np.array(rec.read(lbuf,'d'))
        ibuf = np.array(rec.read(4*lbuf,'B')).reshape(lbuf, 4)
        length = rec.read(1,'i')[0]

        if length < 0: raise StopIteration
        yield buf[:length], ibuf[:length]

def list_integrals(*args, **kwargs):
    """ List two-electron spin-orbit integrals in file """

    for buf, ibuf in list_buffers(*args, **kwargs):
        for g, ig in zip(buf, ibuf):
            yield ig, g


def fockab(Dab, **kwargs):
    """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

    filename = kwargs.get('filename', 'AOTWOINT')
    hfc = kwargs.get('hfc', 1)
    hfx = kwargs.get('hfx', 1)
    f2py = kwargs.get('f2py', True)

    Da, Db = Dab
    D = Da + Db

    if f2py is True:
        try:
            import sirfck
        except(ImportError):
            f2py = False
            print "Warning: non-existent sirfck.so wanted - reverting to python"


    J = matrix(D.shape)
    Ka = matrix(D.shape)
    Kb = matrix(D.shape)

    if f2py:
        for buf, ibuf in list_buffers(filename):
            J, Ka, Kb = sirfck.fckab(
                J, Ka, Kb, Da, Db, buf, ibuf.T
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

    Fa = hfc*J-hfx*Ka
    Fb = hfc*J-hfx*Kb
    return (Fa, Fb)

def fock(D, filename="AOTWOINT", hfc=1, hfx=1, f2py=True):
    """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

    if f2py is True:
        try:
            import sirfck
        except(ImportError):
            f2py = False
            print "Warning: non-existent sirfck.so wanted - reverting to python"

    J = matrix(D.shape)
    K = matrix(D.shape)

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

    print  hfc*J - 0.5*hfx*K
    return hfc*J - 0.5*hfx*K

def semitransform(*args, **kwargs):
    d1, d2 = args
    filename = kwargs.get('file', '/tmp/AOTWOINT')
    m, a = d1.shape
    mmaa = matrix((m, m, a, a))
    for ig, g in list_integrals(filename):
	p, q, r, s = ig
	s, r, q, p = p-1, q-1, r-1, s-1
        if p == q: g *= 0.5
        if r == s: g *= 0.5
	if (p, q) == (r, s): g *= 0.5

	mmaa[:, :, r, s] += d1[:, p].x(d2[:, q])*g
	mmaa[:, :, r, s] += d1[:, q].x(d2[:, p])*g
	mmaa[:, :, s, r] += d1[:, p].x(d2[:, q])*g
	mmaa[:, :, s, r] += d1[:, q].x(d2[:, p])*g
	mmaa[:, :, p, q] += d1[:, r].x(d2[:, s])*g
	mmaa[:, :, q, p] += d1[:, r].x(d2[:, s])*g
	mmaa[:, :, p, q] += d1[:, s].x(d2[:, r])*g
	mmaa[:, :, q, p] += d1[:, s].x(d2[:, r])*g

    return mmaa

if __name__ == "__main__":

    import sys, os, argparse 

    parser = argparse.ArgumentParser()
    parser.add_argument('--file', default='/tmp/AOTWOINT', help='Two-electron integral file')
    parser.add_argument( '--list', action='store_true', help='List integrals on file')

    args = parser.parse_args()

    if args.list:
        print "List integrals"
        for ig, g in list_integrals(args.file):
            print ig, g
