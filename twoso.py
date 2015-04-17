#!/usr/bin/env python
""" This module reads the two-electron spin-orbit file of a Dalton calculation 
    and provides functions for fock matrices

    The structure of the file is
    buf(lbuf), ibuf(lbuf), n
    where lbuf is parameter declared constant, n number integrals in buffer
    For indexi buf[i] is integral and ibuf[i]=(p,q,r,s) is packed 4-index,
    typically one byte per index cartesian components of integralse are 
    intermixed, p=0 marks change of component q is component
"""

import os, sys
import numpy as np
from daltools.util.unformatted import FortranBinary as fb
from daltools.util.full import matrix
import two

def list_integrals(*args, **kwargs):
    """ List two-electron spin-orbit integrals in file """

    for buf, ibuf in two.list_buffers(*args, label='AO2SOINT', **kwargs):
        for g, ig in zip(buf, ibuf):
            if ig[0] == 0:
                comp = ig[1] - 1
            else:
                yield comp, ig, g


def fock(D, component, S_g=1,  **kwargs):
    """ Generate two-electron spin-orbit Fock matrix 
        from integral file AO2SOINT
   """
    
    hfc = kwargs.get('hfc', 1.0)
    hfx = kwargs.get('hfx', 1.0)

#
# These factors for the coulomb part (1 or 2)  are derived from the spin 
# combination  l21 . (s2 + 2s1) which are reversed for singlet and triplet 
# gradients
# The exchange parts are identical
#
    if S_g:
        left = 1
        right = 2
    else:
        left = 2
        right = 1

    nbast, _ = D.shape
    J = matrix((nbast, nbast))
    JL = matrix((nbast, nbast))
    JR = matrix((nbast, nbast))
    K = matrix((nbast, nbast))

    for buf, ibuf in two.list_buffers(label='AO2SOINT', **kwargs):
        for g, ig in zip(buf, ibuf):
            if ig[0] == 0:
                comp = "*xyz"[ig[1]]
            else:
                #print comp, ig, g, component
                if comp != component:
                    continue
                p, q, r, s = ig
                s, r, q, p = (p-1, q-1, r-1, s-1)
                if (p == q): g *= 0.5
                x = 1.5 * g

                JL[r, s] += left * g * (D[p, q] + D[q, p])
                JL[s, r] -= left * g * (D[p, q] + D[q, p])
                JR[p, q] += right* g * (D[r, s] - D[s, r])
                JR[q, p] += right* g * (D[r, s] - D[s, r])

                K[p, s] -= x * D[r, q]
                K[s, p] += x * D[q, r]
                K[p, r] += x * D[s, q]
                K[r, p] -= x * D[q, s]
                K[q, s] -= x * D[r, p]
                K[s, q] += x * D[p, r]
                K[q, r] += x * D[s, p]
                K[r, q] -= x * D[p, s]

    J = JR+JL
    F = hfc*J+hfx*K
    return F

def fockab(Da, Db, component, **kwargs):
    """ Generate two-electron spin-orbit Fock(alpha,beta) matrix from 
        integral file AO2SOINT
        Input: component 'x', 'y', or 'z'
               Density matrices, tuple (Da, Db)
               AO intergral file, FortranBinary object, positioned
        Output Fock matrics, tuple(Fa, Fb)
    """

    hfc = kwargs.get('hfc', 1.0)
    hfx = kwargs.get('hfx', 1.0)
    
    Fa = matrix(Da.shape)
    Fb = matrix(Db.shape)
    

    for buf, ibuf in two.list_buffers(label='AO2SOINT', **kwargs):
        for g, ig in zip(buf, ibuf):
            if ig[0] == 0:
                comp = "*xyz"[ig[1]]
            else:
                #print comp, ig, g, component
                if comp != component:
                    continue
                p, q, r, s = ig
                s, r, q, p = (p-1, q-1, r-1, s-1)
                if (p == q): g *= 0.5
                j = hfc*g
                x = 3*hfx*g

                #Fa[p,q]=(pq|rs)(2D+(r,s) + D-(r,s))
                Fapq = j*(3*(Da[r, s] - Da[s, r]) + Db[r, s] - Db[s, r])
                Fa[p, q] += Fapq
                Fa[q, p] += Fapq

                #Fb[p,q]=(pq|rs)(-2D+(r,s) + D-(r,s))
                Fbpq = j*(-(Da[r, s] - Da[s, r]) - 3*(Db[r, s] - Db[s, r]))
                Fb[p, q] += Fbpq
                Fb[q, p] += Fbpq

                #Fa[r,s]=(pq|rs)(2D-(p,p) + D+(p,q))
                Fars = j*(3*(Da[p, q] + Da[q, p]) - Db[p, q] - Db[q, p])
                Fa[r, s] += Fars
                Fa[s, r] -= Fars

                #Fb[r,s]=(pq|rs)(2D+(p,q) - D-(p,q))
                Fbrs = j*(Da[p, q] + Da[q, p] - 3*(Db[p, q] + Db[q, p]))
                Fb[r, s] += Fbrs
                Fb[s, r] -= Fbrs

                Fa[p, s] -= x*Da[r, q]
                Fa[q, s] -= x*Da[r, p]
                Fa[p, r] += x*Da[s, q]
                Fa[q, r] += x*Da[s, p]

                Fb[p, s] += x*Db[r, q]
                Fb[q, s] += x*Db[r, p]
                Fb[p, r] -= x*Db[s, q]
                Fb[q, r] -= x*Db[s, p]

                Fa[r, q] -= x*Da[p, s]
                Fa[r, p] -= x*Da[q, s]
                Fa[s, q] += x*Da[p, r]
                Fa[s, p] += x*Da[q, r]

                Fb[r, q] += x*Db[p, s]
                Fb[r, p] += x*Db[q, s]
                Fb[s, q] -= x*Db[p, r]
                Fb[s, p] -= x*Db[q, r]

    return (Fa, Fb)


if __name__ == "__main__":
    import optparse 
    parser = optparse.OptionParser()
    parser.add_option(
        '-d','--directory', dest='dir', default='/tmp', 
        help='Directory containing Dalton job files'
        )
    parser.add_option(
        '-f','--fock', dest='fock', default=False, action='store_true', 
        help='Calculate the [inactive] Fock matrix'
        )
    parser.add_option(
        '-g','--grad', dest='grad', default=1, type=int, 
        help='Calculate the [inactive] Fock matrix'
        )
    parser.add_option(
        '-l','--list', dest='list', default=False, action='store_true', 
        help='List integrals on file'
        )
    parser.add_option(
        '-o','--occ', dest='occ', default=None, 
        help='List integrals on file'
        )
    parser.add_option(
        '-c','--comp', dest='comp', default=None, 
        help='List integrals on file'
        )
    parser.add_option(
        '-a','--fock-ab', dest='fock_ab', default=False, action='store_true', 
        help='Calculate the alpha,beta Fock matrices'
        )
    parser.add_option(
        '-J','--hfc', dest='hfc', default=1.0, type=float, 
        help='Hartree-Fock coulomb factor'
        )
    parser.add_option(
        '-K','--hfx', dest='hfx', default=1.0, type=float, 
        help='Hartree-Fock exchange factor'
        )
    parser.add_option(
        '-S','--spin-density', dest='S', default=1, type=int, 
        help='Choose input density (0,1)'
        )

    opt, arg = parser.parse_args(sys.argv[1:])

    #
    # Get ao basis dimension from one-integral file
    #
    ao2soint = os.path.sep.join([opt.dir, "AO2SOINT"])


    if opt.list:
        print "List integrals"
        for c, ig, g in list_integrals(ao2soint, "AO2SOINT"):
            print c, ig, g

    if opt.fock:
        if opt.grad not in  [0, 1]:
            print "Valid opt: 0,1"
            sys.exit(1)
        if opt.grad:
            print "Get triplet Fock matrix"
        else:
            print "Get singlet Fock matrix"

    if opt.occ:
        occ = [int(i) for i in opt.occ.split(",")]
        print occ


    if opt.fock or opt.fock_ab:
#
# MO from from ifc
#
        SIRIUS_RST = os.path.sep.join([opt.dir, "SIRIUS.RST"])
        from sirrst import sirrst
        rst = sirrst(SIRIUS_RST)
        cmo = rst.cmo.unblock()

        if opt.occ:
            na, nb = occ
        else:
            na = nb = cmo.shape[1]

        cmoa = cmo[:, :na]
        cmob = cmo[:, :nb]

        Da = cmoa*cmoa.T
        Db = cmob*cmob.T

    if opt.fock:

        D_S = Da + opt.S*Db
        F_Sg = fock(opt.comp, opt.grad, D_S, 
            filename=ao2soint,  hfc=opt.hfc, hfx=opt.hfx)
        print D_S, F_Sg


    if opt.fock_ab:

        Fa, Fb = fockab(opt.comp, (Da, Db), ao2soint, hfc=opt.hfc, hfx=opt.hfx)
        print  "Fa", Fa, "Fb", Fb


    sys.exit(0)
