C
      SUBROUTINE FCKAB(J,KA,KB,DA,DB,BUF,IBUF,NB,ND,LENGTH)
      IMPLICIT NONE
CF2PY INTENT(IN,OUT) J,KA,KB
CF2PY INTENT(IN) DA,DB,BUF,IBUF
      INTEGER NB,ND,IBUF,LENGTH
      DOUBLE PRECISION J,KA,KB,DA,DB,BUF
      DIMENSION J(NB,NB,ND), KA(NB,NB,ND), KB(NB,NB,ND),
     & DA(NB,NB,ND), DB(NB,NB,ND),
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
      INTEGER P, Q, R, S, INTEGR, I
C
C NOTE: Reals and logicals should appear at the end.
C
      DO 100 INTEGR = 1, LENGTH
         P = IBUF(1,INTEGR)
         Q = IBUF(2,INTEGR)
         R = IBUF(3,INTEGR)
         S = IBUF(4,INTEGR)
         GINT = BUF(INTEGR)
         IF (P.EQ.Q) GINT = GINT / 2
         IF (R.EQ.S) GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C coulomb:
         DO I=1, ND
             FADD = GINT*(DA(R,S,I)+DB(R,S,I)+DA(S,R,I)+DB(S,R,I))
             J(P,Q,I) = J(P,Q,I) + FADD
             J(Q,P,I) = J(Q,P,I) + FADD
             FADD = GINT*(DA(P,Q,I)+DB(P,Q,I)+DA(Q,P,I)+DB(Q,P,I))
             J(R,S,I) = J(R,S,I) + FADD
             J(S,R,I) = J(S,R,I) + FADD
C exchange:
             KA(P,S,I) = KA(P,S,I) + GINT*DA(R,Q,I)
             KA(P,R,I) = KA(P,R,I) + GINT*DA(S,Q,I)
             KA(Q,S,I) = KA(Q,S,I) + GINT*DA(R,P,I)
             KA(Q,R,I) = KA(Q,R,I) + GINT*DA(S,P,I)
             KA(R,Q,I) = KA(R,Q,I) + GINT*DA(P,S,I)
             KA(S,Q,I) = KA(S,Q,I) + GINT*DA(P,R,I)
             KA(R,P,I) = KA(R,P,I) + GINT*DA(Q,S,I)
             KA(S,P,I) = KA(S,P,I) + GINT*DA(Q,R,I)
             KB(P,S,I) = KB(P,S,I) + GINT*DB(R,Q,I)
             KB(P,R,I) = KB(P,R,I) + GINT*DB(S,Q,I)
             KB(Q,S,I) = KB(Q,S,I) + GINT*DB(R,P,I)
             KB(Q,R,I) = KB(Q,R,I) + GINT*DB(S,P,I)
             KB(R,Q,I) = KB(R,Q,I) + GINT*DB(P,S,I)
             KB(S,Q,I) = KB(S,Q,I) + GINT*DB(P,R,I)
             KB(R,P,I) = KB(R,P,I) + GINT*DB(Q,S,I)
             KB(S,P,I) = KB(S,P,I) + GINT*DB(Q,R,I)
        END DO
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
