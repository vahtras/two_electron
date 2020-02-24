/* sirfck.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"


/* Subroutine */ int fckab_(doublereal *j, doublereal *ka, doublereal *kb, 
	doublereal *da, doublereal *db, doublereal *buf, integer *ibuf, 
	integer *nb, integer *nd, integer *length)
{
    /* System generated locals */
    integer j_dim1, j_dim2, j_offset, ka_dim1, ka_dim2, ka_offset, kb_dim1, 
	    kb_dim2, kb_offset, da_dim1, da_dim2, da_offset, db_dim1, db_dim2,
	     db_offset, i__1, i__2;

    /* Local variables */
    static integer i__, p, q, r__, s;
    static doublereal fadd, gint;
    static integer integr;

/* F2PY INTENT(IN,OUT) J,KA,KB */
/* F2PY INTENT(IN) DA,DB,BUF,IBUF */

/* based on FCKD03  - input alpha-beta Fock matrices OV */
/* Henrik Koch and Trygve Helgaker 18-NOV-1991. */
/* 970303-tsaue : index permutation */
/* 941011-hjaaj: renamed from FCKDI1 to FCKD03 */
/* DFT modifications T. Helgaker */

/* FILE: priunit.h */

/* NOTE: Reals and logicals should appear at the end. */

    /* Parameter adjustments */
    db_dim1 = *nb;
    db_dim2 = *nb;
    db_offset = 1 + db_dim1 * (1 + db_dim2);
    db -= db_offset;
    da_dim1 = *nb;
    da_dim2 = *nb;
    da_offset = 1 + da_dim1 * (1 + da_dim2);
    da -= da_offset;
    kb_dim1 = *nb;
    kb_dim2 = *nb;
    kb_offset = 1 + kb_dim1 * (1 + kb_dim2);
    kb -= kb_offset;
    ka_dim1 = *nb;
    ka_dim2 = *nb;
    ka_offset = 1 + ka_dim1 * (1 + ka_dim2);
    ka -= ka_offset;
    j_dim1 = *nb;
    j_dim2 = *nb;
    j_offset = 1 + j_dim1 * (1 + j_dim2);
    j -= j_offset;
    ibuf -= 5;
    --buf;

    /* Function Body */
    i__1 = *length;
    for (integr = 1; integr <= i__1; ++integr) {
	p = ibuf[(integr << 2) + 1];
	q = ibuf[(integr << 2) + 2];
	r__ = ibuf[(integr << 2) + 3];
	s = ibuf[(integr << 2) + 4];
	gint = buf[integr];
	if (p == q) {
	    gint /= 2;
	}
	if (r__ == s) {
	    gint /= 2;
	}
	if (p == r__ && s == q) {
	    gint /= 2;
	}
/* coulomb: */
	i__2 = *nd;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fadd = gint * (da[r__ + (s + i__ * da_dim2) * da_dim1] + db[r__ + 
		    (s + i__ * db_dim2) * db_dim1] + da[s + (r__ + i__ * 
		    da_dim2) * da_dim1] + db[s + (r__ + i__ * db_dim2) * 
		    db_dim1]);
	    j[p + (q + i__ * j_dim2) * j_dim1] += fadd;
	    j[q + (p + i__ * j_dim2) * j_dim1] += fadd;
	    fadd = gint * (da[p + (q + i__ * da_dim2) * da_dim1] + db[p + (q 
		    + i__ * db_dim2) * db_dim1] + da[q + (p + i__ * da_dim2) *
		     da_dim1] + db[q + (p + i__ * db_dim2) * db_dim1]);
	    j[r__ + (s + i__ * j_dim2) * j_dim1] += fadd;
	    j[s + (r__ + i__ * j_dim2) * j_dim1] += fadd;
/* exchange: */
	    ka[p + (s + i__ * ka_dim2) * ka_dim1] += gint * da[r__ + (q + i__ 
		    * da_dim2) * da_dim1];
	    ka[p + (r__ + i__ * ka_dim2) * ka_dim1] += gint * da[s + (q + i__ 
		    * da_dim2) * da_dim1];
	    ka[q + (s + i__ * ka_dim2) * ka_dim1] += gint * da[r__ + (p + i__ 
		    * da_dim2) * da_dim1];
	    ka[q + (r__ + i__ * ka_dim2) * ka_dim1] += gint * da[s + (p + i__ 
		    * da_dim2) * da_dim1];
	    ka[r__ + (q + i__ * ka_dim2) * ka_dim1] += gint * da[p + (s + i__ 
		    * da_dim2) * da_dim1];
	    ka[s + (q + i__ * ka_dim2) * ka_dim1] += gint * da[p + (r__ + i__ 
		    * da_dim2) * da_dim1];
	    ka[r__ + (p + i__ * ka_dim2) * ka_dim1] += gint * da[q + (s + i__ 
		    * da_dim2) * da_dim1];
	    ka[s + (p + i__ * ka_dim2) * ka_dim1] += gint * da[q + (r__ + i__ 
		    * da_dim2) * da_dim1];
	    kb[p + (s + i__ * kb_dim2) * kb_dim1] += gint * db[r__ + (q + i__ 
		    * db_dim2) * db_dim1];
	    kb[p + (r__ + i__ * kb_dim2) * kb_dim1] += gint * db[s + (q + i__ 
		    * db_dim2) * db_dim1];
	    kb[q + (s + i__ * kb_dim2) * kb_dim1] += gint * db[r__ + (p + i__ 
		    * db_dim2) * db_dim1];
	    kb[q + (r__ + i__ * kb_dim2) * kb_dim1] += gint * db[s + (p + i__ 
		    * db_dim2) * db_dim1];
	    kb[r__ + (q + i__ * kb_dim2) * kb_dim1] += gint * db[p + (s + i__ 
		    * db_dim2) * db_dim1];
	    kb[s + (q + i__ * kb_dim2) * kb_dim1] += gint * db[p + (r__ + i__ 
		    * db_dim2) * db_dim1];
	    kb[r__ + (p + i__ * kb_dim2) * kb_dim1] += gint * db[q + (s + i__ 
		    * db_dim2) * db_dim1];
	    kb[s + (p + i__ * kb_dim2) * kb_dim1] += gint * db[q + (r__ + i__ 
		    * db_dim2) * db_dim1];
	}
/* L100: */
    }
    return 0;
} /* fckab_ */

/* Subroutine */ int fck_(doublereal *j, doublereal *k, doublereal *dj, 
	doublereal *dk, doublereal *buf, integer *ibuf, integer *n, integer *
	length)
{
    /* System generated locals */
    integer j_dim1, j_offset, k_dim1, k_offset, dj_dim1, dj_offset, dk_dim1, 
	    dk_offset, i__1;

    /* Local variables */
    static integer p, q, r__, s, int__;
    static doublereal fadd, gint;

/* F2PY INTENT(IN,OUT) J,K */
/* F2PY INTENT(IN) DJ,DK,BUF,IBUF */

/* based on FCKD03  - input alpha-beta Fock matrices OV */
/* Henrik Koch and Trygve Helgaker 18-NOV-1991. */
/* 970303-tsaue : index permutation */
/* 941011-hjaaj: renamed from FCKDI1 to FCKD03 */
/* DFT modifications T. Helgaker */

/* FILE: priunit.h */

/* NOTE: Reals and logicals should appear at the end. */

    /* Parameter adjustments */
    dk_dim1 = *n;
    dk_offset = 1 + dk_dim1;
    dk -= dk_offset;
    dj_dim1 = *n;
    dj_offset = 1 + dj_dim1;
    dj -= dj_offset;
    k_dim1 = *n;
    k_offset = 1 + k_dim1;
    k -= k_offset;
    j_dim1 = *n;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    ibuf -= 5;
    --buf;

    /* Function Body */
    i__1 = *length;
    for (int__ = 1; int__ <= i__1; ++int__) {
	p = ibuf[(int__ << 2) + 1];
	q = ibuf[(int__ << 2) + 2];
	r__ = ibuf[(int__ << 2) + 3];
	s = ibuf[(int__ << 2) + 4];
	gint = buf[int__];
	if (p == q) {
	    gint /= 2;
	}
	if (r__ == s) {
	    gint /= 2;
	}
	if (p == r__ && s == q) {
	    gint /= 2;
	}
/* coulomb: */
	fadd = gint * (dj[r__ + s * dj_dim1] + dj[s + r__ * dj_dim1]);
	j[p + q * j_dim1] += fadd;
	j[q + p * j_dim1] += fadd;
	fadd = gint * (dj[p + q * dj_dim1] + dj[q + p * dj_dim1]);
	j[r__ + s * j_dim1] += fadd;
	j[s + r__ * j_dim1] += fadd;
/* exchange: */
/*        GINT = GINT/2 */
	k[p + s * k_dim1] += gint * dk[r__ + q * dk_dim1];
	k[p + r__ * k_dim1] += gint * dk[s + q * dk_dim1];
	k[q + s * k_dim1] += gint * dk[r__ + p * dk_dim1];
	k[q + r__ * k_dim1] += gint * dk[s + p * dk_dim1];
	k[r__ + q * k_dim1] += gint * dk[p + s * dk_dim1];
	k[s + q * k_dim1] += gint * dk[p + r__ * dk_dim1];
	k[r__ + p * k_dim1] += gint * dk[q + s * dk_dim1];
	k[s + p * k_dim1] += gint * dk[q + r__ * dk_dim1];
/* L100: */
    }
    return 0;
} /* fck_ */

