      FUNCTION RAN2(ISEED)!Declan: Changed RAN to RAN2 for gfortran compatability
      IMPLICIT NONE
      INTEGER ISEED,NEVHEP,NRN(2)
      DOUBLE PRECISION RAN2,RAN1,DUMMY,HWRSET,HWRGET
      DATA NEVHEP/0/
      NEVHEP=NEVHEP+1
      NRN(1)=ISEED
      NRN(2)=ABS(ISEED-111111111)
      IF(NEVHEP.EQ.1)DUMMY = HWRSET(NRN)
      RAN2=RAN1(0)
      RETURN
      END

      FUNCTION RAN1(I)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION RAN1,HWRSET,HWRGET
      INTEGER I,ISEED(2),K,IZ,JSEED(2)
      SAVE ISEED
      DATA ISEED/12345,67890/
      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      RAN1=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      RETURN
C-----------------------------------------------------------------------
      ENTRY HWRSET(JSEED)
C-----------------------------------------------------------------------
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
C-----------------------------------------------------------------------
      ENTRY HWRGET(JSEED)
C-----------------------------------------------------------------------
      JSEED(1)=ISEED(1)
      JSEED(2)=ISEED(2)
      RETURN
      END
