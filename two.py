#!/usr/bin/env python
"""Module for getting stuff from Dalton two-electron integral file AOTWOINT"""

import numpy as np
from util.unformatted import FortranBinary as FB
from util.full import matrix

class TwoInt(object):
    def __init__(self, aotwoint):
        self.aotwoint = aotwoint

    def info(self):
        """Extract info block BASINFO on AOTWOINT"""
        _aotwoint = FB(self.aotwoint)
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

    def fock(self, D, hfc=1, hfx=1, f2py=True):
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
            for buf, ibuf in self.list_buffers():
                J, K = sirfck.fck(
                    J, K,  D, D, buf, ibuf.T
                    )
        else:
            for ig, g in self.list_integrals():
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

    def list_integrals(self, *args, **kwargs):
        """ List two-electron spin-orbit integrals in file """

        for buf, ibuf in self.list_buffers(*args, **kwargs):
            for g, ig in zip(buf, ibuf):
                yield ig, g

    def list_buffers(self, label="BASTWOEL"):
        """ Return integral buffers in AOTWOINT"""
        _aofile = FB(self.aotwoint, label=label)

        for rec in _aofile:
            lbuf = (_aofile.reclen-4)/12

            buf = np.array(rec.read(lbuf,'d'))
            ibuf = np.array(rec.read(4*lbuf,'B')).reshape(lbuf, 4)
            length = rec.read(1,'i')[0]

            if length < 0: raise StopIteration
            yield buf[:length], ibuf[:length]


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

    return hfc*J - 0.5*hfx*K

def vb_transform(dens, delta, **kwargs):
    filename = kwargs.get('file', '/tmp/AOTWOINT')
    a, m = dens[0].shape
    H_uvmn = matrix((a, a, m, m))
    for ig, g in list_integrals(filename):
        p, q, r, s = ig
        s, r, q, p = p-1, q-1, r-1, s-1
        if p == q: g *= 0.5
        if r == s: g *= 0.5
        if (p, q) == (r, s): g *= 0.5

        for d1, D1 in zip(dens, delta):
            for d2, D2 in zip(dens, delta):
                H_uvmn += D1[:, q].x(D2[:, s].x(d1[p, :].x(d2[r, :]*g)))
                H_uvmn += D1[:, p].x(D2[:, s].x(d1[q, :].x(d2[r, :]*g)))
                H_uvmn += D1[:, q].x(D2[:, r].x(d1[p, :].x(d2[s, :]*g)))
                H_uvmn += D1[:, p].x(D2[:, r].x(d1[q, :].x(d2[s, :]*g)))
                H_uvmn += D1[:, s].x(D2[:, q].x(d1[r, :].x(d2[p, :]*g)))
                H_uvmn += D1[:, r].x(D2[:, q].x(d1[s, :].x(d2[p, :]*g)))
                H_uvmn += D1[:, s].x(D2[:, p].x(d1[r, :].x(d2[q, :]*g)))
                H_uvmn += D1[:, r].x(D2[:, p].x(d1[s, :].x(d2[q, :]*g)))

            H_uvmn -= D1[:, s].x(D1[:, q].x(d1[p, :].x(d1[r, :]*g)))
            H_uvmn -= D1[:, s].x(D1[:, p].x(d1[q, :].x(d1[r, :]*g)))
            H_uvmn -= D1[:, r].x(D1[:, q].x(d1[p, :].x(d1[s, :]*g)))
            H_uvmn -= D1[:, r].x(D1[:, p].x(d1[q, :].x(d1[s, :]*g)))
            H_uvmn -= D1[:, q].x(D1[:, s].x(d1[r, :].x(d1[p, :]*g)))
            H_uvmn -= D1[:, q].x(D1[:, r].x(d1[s, :].x(d1[p, :]*g)))
            H_uvmn -= D1[:, p].x(D1[:, s].x(d1[r, :].x(d1[q, :]*g)))
            H_uvmn -= D1[:, p].x(D1[:, r].x(d1[s, :].x(d1[q, :]*g)))
        
    return H_uvmn.transpose(0, 2, 1, 3)

def vb_transform2(Dma, Dam, Delta1, Delta2, **kwargs):
    filename = kwargs.get('file', '/tmp/AOTWOINT')
    H_umvn = matrix(Dam[0].shape + Dam[0].shape)
    for ig, g in list_integrals(filename):
        p, q, r, s = ig
        s, r, q, p = p-1, q-1, r-1, s-1
        if p == q: g *= 0.5
        if r == s: g *= 0.5
        if (p, q) == (r, s): g *= 0.5

        for d1a, d2a, D1a, D2a in zip(Dma, Dam, Delta1, Delta2):
            for d1b, d2b, D1b, D2b in zip(Dma, Dam, Delta1, Delta2):
               H_umvn += D1a[p,:].x(d1a[:, q].x(D2b[:, s].x(d2b[r, :]*g)))
               H_umvn += D1a[q,:].x(d1a[:, p].x(D2b[:, s].x(d2b[r, :]*g)))
               H_umvn += D1a[p,:].x(d1a[:, q].x(D2b[:, r].x(d2b[s, :]*g)))
               H_umvn += D1a[q,:].x(d1a[:, p].x(D2b[:, r].x(d2b[s, :]*g)))

               H_umvn += D1a[r,:].x(d1a[:, s].x(D2b[:, q].x(d2b[p, :]*g)))
               H_umvn += D1a[r,:].x(d1a[:, s].x(D2b[:, p].x(d2b[q, :]*g)))
               H_umvn += D1a[s,:].x(d1a[:, r].x(D2b[:, q].x(d2b[p, :]*g)))
               H_umvn += D1a[s,:].x(d1a[:, r].x(D2b[:, p].x(d2b[q, :]*g)))

            H_umvn -= D1a[p, :].x(d1a[:, s].x(D2a[:, q].x(d2a[r, :]*g)))
            H_umvn -= D1a[q, :].x(d1a[:, s].x(D2a[:, p].x(d2a[r, :]*g)))
            H_umvn -= D1a[p, :].x(d1a[:, r].x(D2a[:, q].x(d2a[s, :]*g)))
            H_umvn -= D1a[q, :].x(d1a[:, r].x(D2a[:, p].x(d2a[s, :]*g)))

            H_umvn -= D1a[r, :].x(d1a[:, q].x(D2a[:, s].x(d2a[p, :]*g)))
            H_umvn -= D1a[r, :].x(d1a[:, p].x(D2a[:, s].x(d2a[q, :]*g)))
            H_umvn -= D1a[s, :].x(d1a[:, q].x(D2a[:, r].x(d2a[p, :]*g)))
            H_umvn -= D1a[s, :].x(d1a[:, p].x(D2a[:, r].x(d2a[q, :]*g)))

    return H_umvn




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
