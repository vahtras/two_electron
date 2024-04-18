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

import sqlite3
import sys
from util.full import matrix
from . import core as two
from .eri import Reader

class SpinOrbitReader(Reader):

    def list_integrals(self):
        """List two-electron spin-orbit integrals in file"""

        for buf, ibuf in super().list_buffers(label="AO2SOINT"):
            for g, ig in zip(buf, ibuf):
                if ig[0] == 0:
                    comp = ig[1] - 1
                else:
                    yield comp, ig, g

    def to_sql(self, con):
        """Create SQL table for two-electron spin-orbit integrals"""

        cur = con.cursor()
        cur.execute(
            """CREATE TABLE IF NOT EXISTS ao2soint (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            comp INTEGER,
            p INTEGER,
            q INTEGER,
            r INTEGER,
            s INTEGER,
            value REAL)"""
        )

        # cur.execute("CREATE INDEX IF NOT EXISTS ao2soint_idx ON ao2soint(comp,p,q,r,s)")

        cur.execute("BEGIN TRANSACTION")
        for comp, ig, g in self.list_integrals():
            comp = int(comp)
            p, q, r, s = (int(_) for _ in ig)
            g = float(g)
            record = (comp, p, q, r, s, g)
            print(record)
            cur.execute(
                "INSERT INTO ao2soint(comp,p,q,r,s,value) VALUES(?,?,?,?,?,?)",
                record,
            )
        cur.execute("COMMIT")

    def fock(self, D, component, **kwargs):
        """Generate two-electron spin-orbit Fock matrix
        from integral file AO2SOINT
        """

        hfc = kwargs.get("hfc", 1.0)
        hfx = kwargs.get("hfx", 1.0)

        nbast, _ = D.shape
        J = matrix((nbast, nbast))
        JL = matrix((nbast, nbast))
        JR = matrix((nbast, nbast))
        K = matrix((nbast, nbast))

        for buf, ibuf in two.list_buffers(label="AO2SOINT", **kwargs):
            for g, ig in zip(buf, ibuf):
                if ig[0] == 0:
                    comp = "*xyz"[ig[1]]
                else:
                    # print comp, ig, g, component
                    if comp != component:
                        continue
                    p, q, r, s = ig
                    s, r, q, p = (p - 1, q - 1, r - 1, s - 1)
                    if p == q:
                        g *= 0.5
                    x = 1.5 * g

                    JL[r, s] += g * (D[p, q] + D[q, p])
                    JL[s, r] -= g * (D[p, q] + D[q, p])
                    JR[p, q] += 2 * g * (D[r, s] - D[s, r])
                    JR[q, p] += 2 * g * (D[r, s] - D[s, r])

                    K[p, s] -= x * D[r, q]
                    K[s, p] += x * D[q, r]
                    K[p, r] += x * D[s, q]
                    K[r, p] -= x * D[q, s]
                    K[q, s] -= x * D[r, p]
                    K[s, q] += x * D[p, r]
                    K[q, r] += x * D[s, p]
                    K[r, q] -= x * D[p, s]

        J = JR + JL
        F = hfc * J + hfx * K
        return F


    def fockab(self, Da, Db, component, **kwargs):
        """Generate two-electron spin-orbit Fock(alpha,beta) matrix from
        integral file AO2SOINT
        Input: component 'x', 'y', or 'z'
               Density matrices, tuple (Da, Db)
               AO integral file, FortranBinary object, positioned
        Output Fock matrics, tuple(Fa, Fb)
        """

        hfc = kwargs.get("hfc", 1.0)
        hfx = kwargs.get("hfx", 1.0)

        Fa = matrix(Da.shape)
        Fb = matrix(Db.shape)

        for buf, ibuf in super().list_buffers(label="AO2SOINT", **kwargs):
            for g, ig in zip(buf, ibuf):
                if ig[0] == 0:
                    comp = "*xyz"[ig[1]]
                else:
                    # print comp, ig, g, component
                    if comp != component:
                        continue
                    p, q, r, s = ig
                    s, r, q, p = (p - 1, q - 1, r - 1, s - 1)
                    if p == q:
                        g *= 0.5
                    j = hfc * g
                    x = 3 * hfx * g

                    # Fa[p,q]=(pq|rs)(2D+(r,s) + D-(r,s))
                    Fapq = j * (3 * (Da[r, s] - Da[s, r]) + Db[r, s] - Db[s, r])
                    Fa[p, q] += Fapq
                    Fa[q, p] += Fapq

                    # Fb[p,q]=(pq|rs)(-2D+(r,s) + D-(r,s))
                    Fbpq = j * (-(Da[r, s] - Da[s, r]) - 3 * (Db[r, s] - Db[s, r]))
                    Fb[p, q] += Fbpq
                    Fb[q, p] += Fbpq

                    # Fa[r,s]=(pq|rs)(2D-(p,p) + D+(p,q))
                    Fars = j * (3 * (Da[p, q] + Da[q, p]) - Db[p, q] - Db[q, p])
                    Fa[r, s] += Fars
                    Fa[s, r] -= Fars

                    # Fb[r,s]=(pq|rs)(2D+(p,q) - D-(p,q))
                    Fbrs = j * (Da[p, q] + Da[q, p] - 3 * (Db[p, q] + Db[q, p]))
                    Fb[r, s] += Fbrs
                    Fb[s, r] -= Fbrs

                    Fa[p, s] -= x * Da[r, q]
                    Fa[q, s] -= x * Da[r, p]
                    Fa[p, r] += x * Da[s, q]
                    Fa[q, r] += x * Da[s, p]

                    Fb[p, s] += x * Db[r, q]
                    Fb[q, s] += x * Db[r, p]
                    Fb[p, r] -= x * Db[s, q]
                    Fb[q, r] -= x * Db[s, p]

                    Fa[r, q] -= x * Da[p, s]
                    Fa[r, p] -= x * Da[q, s]
                    Fa[s, q] += x * Da[p, r]
                    Fa[s, p] += x * Da[q, r]

                    Fb[r, q] += x * Db[p, s]
                    Fb[r, p] += x * Db[q, s]
                    Fb[s, q] -= x * Db[p, r]
                    Fb[s, p] -= x * Db[q, r]

        return (Fa, Fb)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ao2soint", default="AO2SOINT", help="Integral file"
    )
    parser.add_argument(
        "-l", "--list", action="store_true", help="List integrals on file"
    )
    parser.add_argument('--to-sql', action="store_true", help="Create SQL table")

    args = parser.parse_args()

    #
    # Get ao basis dimension from one-integral file
    #

    reader = SpinOrbitReader(args.ao2soint)
    if args.to_sql:
        con = sqlite3.connect('ao2soint.db')
        reader.to_sql(con)
    if args.list:
        print("List integrals")
        for c, ig, g in reader.list_integrals():
            print(c, ig, g)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
