C
      SUBROUTINE FCKAB(J,KA,KB,DA,DB,BUF,IBUF,N,LENGTH)
      IMPLICIT NONE
CF2PY INTENT(IN,OUT) J,KA,KB
CF2PY INTENT(IN) DA,DB,BUF,IBUF
      INTEGER N,IBUF,LENGTH
      DOUBLE PRECISION J,KA,KB,DA,DB,BUF
      DIMENSION J(N,N),KA(N,N),KB(N,N),DA(N,N),DB(N,N),
     & BUF(LENGTH),IBUF(4,LENGTH)
C
C based on FCKD03  - input alpha-beta Fock matrices OV
C Henrik Koch and Trygve Helgaker 18-NOV-1991.
C 970303-tsaue : index permutation
C 941011-hjaaj: renamed from FCKDI1 to FCKD03
C DFT modifications T. Helgaker
C
C FILE: priunit.h
      DOUBLE PRECISION GINT,FADD
      INTEGER P, Q, R, S, INT
C
C NOTE: Reals and logicals should appear at the end.
C
      DO 100 INT = 1, LENGTH
         P = IBUF(1,INT)
         Q = IBUF(2,INT)
         R = IBUF(3,INT)
         S = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q) GINT = GINT / 2
         IF (R.EQ.S) GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C coulomb:
         FADD = GINT*(DA(R,S)+DB(R,S)+DA(S,R)+DB(S,R))
         J(P,Q) = J(P,Q) + FADD
         J(Q,P) = J(Q,P) + FADD
         FADD = GINT*(DA(P,Q)+DB(P,Q)+DA(Q,P)+DB(Q,P))
         J(R,S) = J(R,S) + FADD
         J(S,R) = J(S,R) + FADD
C exchange:
         KA(P,S) = KA(P,S) + GINT*DA(R,Q)
         KA(P,R) = KA(P,R) + GINT*DA(S,Q)
         KA(Q,S) = KA(Q,S) + GINT*DA(R,P)
         KA(Q,R) = KA(Q,R) + GINT*DA(S,P)
         KA(R,Q) = KA(R,Q) + GINT*DA(P,S)
         KA(S,Q) = KA(S,Q) + GINT*DA(P,R)
         KA(R,P) = KA(R,P) + GINT*DA(Q,S)
         KA(S,P) = KA(S,P) + GINT*DA(Q,R)
         KB(P,S) = KB(P,S) + GINT*DB(R,Q)
         KB(P,R) = KB(P,R) + GINT*DB(S,Q)
         KB(Q,S) = KB(Q,S) + GINT*DB(R,P)
         KB(Q,R) = KB(Q,R) + GINT*DB(S,P)
         KB(R,Q) = KB(R,Q) + GINT*DB(P,S)
         KB(S,Q) = KB(S,Q) + GINT*DB(P,R)
         KB(R,P) = KB(R,P) + GINT*DB(Q,S)
         KB(S,P) = KB(S,P) + GINT*DB(Q,R)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE FCK(J,K,DJ,DK,BUF,IBUF,N,LENGTH)
      IMPLICIT NONE
CF2PY INTENT(IN,OUT) J,K
CF2PY INTENT(IN) DJ,DK,BUF,IBUF
      INTEGER N,IBUF,LENGTH
      DOUBLE PRECISION J,K,DJ,DK,BUF
      DIMENSION J(N,N),K(N,N),DJ(N,N),DK(N,N),
     & BUF(LENGTH),IBUF(4,LENGTH)
C
C based on FCKD03  - input alpha-beta Fock matrices OV
C Henrik Koch and Trygve Helgaker 18-NOV-1991.
C 970303-tsaue : index permutation
C 941011-hjaaj: renamed from FCKDI1 to FCKD03
C DFT modifications T. Helgaker
C
C FILE: priunit.h
      DOUBLE PRECISION GINT,FADD
      INTEGER P, Q, R, S, INT
C
C NOTE: Reals and logicals should appear at the end.
C
      DO 100 INT = 1, LENGTH
         P = IBUF(1,INT)
         Q = IBUF(2,INT)
         R = IBUF(3,INT)
         S = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q) GINT = GINT / 2
         IF (R.EQ.S) GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C coulomb:
         FADD = GINT*(DJ(R,S)+DJ(S,R))
         J(P,Q) = J(P,Q) + FADD
         J(Q,P) = J(Q,P) + FADD
         FADD = GINT*(DJ(P,Q)+DJ(Q,P))
         J(R,S) = J(R,S) + FADD
         J(S,R) = J(S,R) + FADD
C exchange:
C        GINT = GINT/2
         K(P,S) = K(P,S) + GINT*DK(R,Q)
         K(P,R) = K(P,R) + GINT*DK(S,Q)
         K(Q,S) = K(Q,S) + GINT*DK(R,P)
         K(Q,R) = K(Q,R) + GINT*DK(S,P)
         K(R,Q) = K(R,Q) + GINT*DK(P,S)
         K(S,Q) = K(S,Q) + GINT*DK(P,R)
         K(R,P) = K(R,P) + GINT*DK(Q,S)
         K(S,P) = K(S,P) + GINT*DK(Q,R)
  100 CONTINUE
      RETURN
      END
