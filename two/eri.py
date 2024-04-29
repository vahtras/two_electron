import pathlib
import sys
import sqlite3

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
            integrals = zip(buf, ibuf)
            for g, ig in integrals:
                yield tuple(ig), g

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

    def fock(self, D, hfc=1, hfx=1):
        """
        Routine for building J and K Fock matrices reading AOTWOINT
        """

        J = np.zeros(D.shape)
        K = np.zeros(D.shape)

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

        Fab = self.fock_builder(Dab, **kwargs)
        return Fab

    def fock_builder(self, Dab, **kwargs):

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




class FReader(Reader):

    def fock(self, D, hfc=1, hfx=1):
        """Routine for building alpha and beta fock matrix by reading AOTWOINT"""

        J = np.zeros(D.shape)
        K = np.zeros(D.shape)

        for buf, ibuf in self.list_buffers():
            J, K = sirfck.fck(
                J, K,  D, D, buf, ibuf.T
                )

        return hfc*J - 0.5*hfx*K

    def fock_builder(self, Dab, **kwargs):

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
        for buf, ibuf in self.list_buffers():
            Js, Kas, Kbs = sirfck.fckab(
                Js, Kas, Kbs, Das, Dbs, buf, ibuf.T
                )

        Fas = (hfc*Js[:, :, i] - hfx*Kas[:, :, i] for i in range(len(Dab)))
        Fbs = (hfc*Js[:, :, i] - hfx*Kbs[:, :, i] for i in range(len(Dab)))

        Fab = tuple((Fa, Fb) for Fa, Fb in zip(Fas, Fbs))
        return Fab

class SQLReader(Reader):
    def __init__(self, filename, db=None):
        super().__init__(filename)
        if db:
            self.db = db
        else:
            self.db = pathlib.Path(self.filename).with_suffix(".db")
        self.con = sqlite3.connect(self.db)

    def insert_integrals(self):
        self.to_sql(self.con)

    def to_sql(self, con):
        """Create SQL table for two-electron spin-orbit integrals"""

        cur = con.cursor()
        cur.execute(
            """
            DROP TABLE IF EXISTS aotwoint;
            """
        )
        cur.execute(
            """
            CREATE TABLE aotwoint (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            p INTEGER,
            q INTEGER,
            r INTEGER,
            s INTEGER,
            value REAL);
            """
        )

        # cur.execute("CREATE INDEX IF NOT EXISTS aotwoint_idx ON aotwoint(comp,p,q,r,s)")

        cur.execute("BEGIN TRANSACTION")
        for ig, g in super().list_integrals():
            p, q, r, s = (int(_) for _ in ig)
            g = float(g)
            s, r, q, p = (p-1, q-1, r-1, s-1)
            record = (p, q, r, s, g)
            cur.execute(
                "INSERT INTO aotwoint(p,q,r,s,value) VALUES(?,?,?,?,?)",
                record,
            )
        cur.execute("COMMIT")

    def list_integrals(self):
        """
        List two-electron spin-orbit integrals in file
        """
        cur = self.con.cursor()
        cur.execute("SELECT p, q, r, s, value FROM aotwoint")
        for ig in cur:
            yield ig[:4], ig[4]

    def insert_density(self, D):
        cur = self.con.cursor()
        cur.execute("DROP TABLE IF EXISTS density;")
        cur.execute("CREATE TABLE density (t INTEGER, u INTEGER, value REAL);")
        for (t, u), value in np.ndenumerate(D):
            if abs(value) < 1e-10:
                continue
            cur.execute("INSERT INTO density VALUES(?, ?, ?);", (t, u, value))
        cur.execute("COMMIT")

    def list_density(self):
        cur = self.con.cursor()
        cur.execute("SELECT t, u, value FROM density")
        for pq in cur:
            yield pq[:2], pq[2]

    def fock(self, D, hfc=1, hfx=1):
        self.insert_density(D)
        cur = self.con.cursor()
        #
        # F(p,q) = (pq|rs)D(r,s)
        #
        if True:
            records1a = [
                (rec[1],rec[2],rec[9])
                for rec in cur.execute(
                    """
                    SELECT *,SUM(aotwoint.value*density.value) AS fock
                    FROM aotwoint JOIN density
                    ON (r=t AND s=u)
                    GROUP BY p, q;
                    """
                )
            ] 
            records1b = [
                (rec[1],rec[2],rec[9])
                for rec in cur.execute(
                    """
                    SELECT *,SUM(aotwoint.value*density.value) AS fock
                    FROM aotwoint JOIN density
                    ON (r=u AND s=t)
                    WHERE r != s
                    GROUP BY p, q;
                    """
                )
            ] 
            f1 = records_to_array(records1a, D.shape, symmetrize=True)
            f1 += records_to_array(records1b, D.shape, symmetrize=True)
        else:
            records1 = [
                (rec[1],rec[2],rec[9])
                for rec in cur.execute(
                    """
                    SELECT *,SUM(aotwoint.value*density.value) AS fock
                    FROM aotwoint JOIN density
                    ON (r=t AND s=u) OR (r=u AND s=t)
                    GROUP BY p, q;
                    """
                )
            ] 
            f1 = records_to_array(records1, D.shape, symmetrize=True)
        #
        # F(r,s) = (pq|rs)D(p,q)
        # skipping diagonal part (pq|pq) to avoid double counting
        #
        if True:
            records2a = [
                (rec[3],rec[4],rec[9])
                for rec in cur.execute(
                    """
                    SELECT *,SUM(aotwoint.value*density.value) AS fock
                    FROM aotwoint JOIN density
                    ON (p=t AND q=u)
                    WHERE NOT (p = r AND q = s)
                    GROUP BY r, s;
                    """
                )
            ]
            records2b = [
                (rec[3],rec[4],rec[9])
                for rec in cur.execute(
                    """
                    SELECT *,SUM(aotwoint.value*density.value) AS fock
                    FROM aotwoint JOIN density
                    ON (p=u AND q=t)
                    WHERE NOT (p = r AND q = s) AND p != q
                    GROUP BY r, s;
                    """
                )
            ]
            f2 = records_to_array(records2a, D.shape, symmetrize=True)
            f2 += records_to_array(records2b, D.shape, symmetrize=True)
        else:
            records2 = [
                (rec[3],rec[4],rec[9])
                for rec in cur.execute(
                    """
                    SELECT *,SUM(aotwoint.value*density.value) AS fock
                    FROM aotwoint JOIN density
                    ON (p=t AND q=u) OR (p=u AND q=t)
                    WHERE (p != r OR q != s)
                    GROUP BY r, s;
                    """
                )
            ]
            f2 = records_to_array(records2, D.shape, symmetrize=True)
        j = f1 + f2


        # exchange 
        # K(p, s) = (pq|rs)D(r, q)
        records_ps = [
            (rec[1],rec[4],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (r=t AND q=u)
                GROUP BY p, s;
                """
            )
        ]
        k = records_to_array(records_ps, D.shape)

        # K(p, r) = (pq|sr)D(s, q)
        records_pr = [
            (rec[1],rec[3],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (s=t AND q=u)
                WHERE r != s
                GROUP BY p, r;
                """
            )
        ]
        k += records_to_array(records_pr, D.shape)

        # K(q, s) = (qp|rs)D(r, p)
        records_qs = [
            (rec[2],rec[4],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (r=t AND p=u)
                where p != q
                GROUP BY q, s;
                """
            )
        ]
        k += records_to_array(records_qs, D.shape)

        # K(q, r) = (qp|sr)D(s, p)
        records_qr = [
            (rec[2],rec[3],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (s=t AND p=u)
                WHERE p != q AND r != s
                GROUP BY q, r;
                """
            ) ]
        k += records_to_array(records_qr, D.shape)

        #K(r,q) = (pq|rs)D(p,s)
        records_rq = [
            (rec[3],rec[2],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (p=t AND s=u)
                WHERE NOT (p = r AND q = s)
                GROUP BY r, q;
                """
            )
        ]
        k += records_to_array(records_rq, D.shape)

        #K(s,q) = (pq|sr)D(p,r)
        records_sq = [
            (rec[4],rec[2],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (p=t AND r=u)
                WHERE NOT(p = r AND q = s) AND (r != s)
                GROUP BY s, q;
                """
            )
        ]
        k += records_to_array(records_sq, D.shape)

        #K(r,p) = (qp|rs)D(q,s)
        records_rp = [
            (rec[3],rec[1],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (q=t AND s=u)
                WHERE NOT(p = r AND q = s) AND (p != q)
                GROUP BY r, p;
                """
            )
        ]
        k += records_to_array(records_rp, D.shape)

        #K(s,p) = (qp|sr)D(q,r)
        records_sp = [
            (rec[4],rec[1],rec[9])
            for rec in cur.execute(
                """
                SELECT *,SUM(aotwoint.value*density.value) AS fock 
                FROM aotwoint JOIN density 
                ON (q=t AND r=u)
                WHERE NOT(p = r AND q = s) AND (p != q) AND (r != s)
                GROUP BY s, p;
                """
            )
        ]
        k += records_to_array(records_sp, D.shape)

        f = hfc*j - 0.5*hfx*k
        return f


def records_to_array(records, shape, symmetrize=False):
    f = np.zeros(shape)
    for p, q, value in records:
        f[p, q] = value
        if symmetrize:
            f[q, p] = value
    return f


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
    parser.add_argument('--to-sql', default=None, help="Create SQL table")

    args = parser.parse_args()

    reader = Reader(args.file)
    if args.list:
        print("List integrals")
        for ig, g in reader.list_integrals():
            print(ig, g)

    if args.to_sql:
        SQLReader(args.file, db=args.to_sql).insert_integrals()

if __name__ == "__main__":
    sys.exit(main())
