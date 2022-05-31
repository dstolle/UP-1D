      MODULE GCONTROL
       INTEGER, PARAMETER :: NNG = 1000, NELT  = 1000
       INTEGER, PARAMETER :: NFD = 2000, MSTIF = 100000 
       INTEGER :: NEL,NNOD,NVAR,NNODEL,NVEL,NNET,LBAND,NMAT,NVA,IFix
	 DOUBLE PRECISION :: Rmin,Rmax,DL,TR 
       CHARACTER (LEN = 5) :: INFILE
      END MODULE GCONTROL

	MODULE MPROPERTY
      DOUBLE PRECISION :: GMOD,EMOD,ANV,GRR,CF,XKAP
	END MODULE MPROPERTY

      MODULE MPOINT
        INTEGER, PARAMETER :: IPNT=1000
        DOUBLE PRECISION, DIMENSION (IPNT)::GAMA,FBAR,UW,PT,RC,RP
        DOUBLE PRECISION, DIMENSION (IPNT)::SXX,SYY,SXY,SZZ,UR
      END MODULE MPOINT

C*****************************************************************

      PROGRAM UP1D
      USE GCONTROL; USE MPROPERTY; USE MPOINT			
C       ************************************************************
C       *                                                          *
C       *      2-NODED STRESS ANALYSIS FINITE ELEMENT PROGRAM      *
C       *                 WITH P CONSTANT, SIJ BILINEAR            *
C       *                    This is linear elastic                *
C       *                         (March 2022)                     *
C       ************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION (NFD) :: DISP, TLOAD
      DOUBLE PRECISION :: GSTIF(MSTIF),RG(NNG),DET
      INTEGER :: NCN,ICO(3,NELT),IX(NFD),IP(NFD)

      WRITE(*,'(A,\)') ' INPUT FILE NAME = '
      READ(*,'(A)') INFILE

      OPEN(2,FILE=INFILE//'.DIS',FORM='FORMATTED')
      OPEN(4,FILE=INFILE//'.STR',FORM='FORMATTED')   ! Stresses for consolidation
      OPEN(5,FILE=INFILE//'.DAT',FORM='FORMATTED')
      OPEN(6,FILE=INFILE//'.OUT',FORM='FORMATTED')

C** INPUT THE INFORMATION NEEDED FOR THE REQUESTED FINITE ELEMENT ANALYSIS

      CALL INPUT(RG,ICO,IX)
	
C** INITIALIZES MOST ARRAYS INCLUDING INITIAL STRESSES

      TLOAD(1:NNET) = 0.D0;  TLOAD(1) = TR*Rmin;
      DISP(1:NNET) = TLOAD(1:NNET);
      SXX(1:NEL) = 0.D0;  SYY(1:NEL) = 0.D0;     
      SXY(1:NEL) = 0.D0;  SZZ(1:NEL) = 0.D0;
      UW(1:NEL) = 0.D0;	     

      DET=1.D-8
      CALL MATRX(GSTIF,ICO,RG,IX)

      CALL GBAND(GSTIF,DISP,1,IP,DET,NCN,LBAND,NNET)
      IF(DET.LE.0.D0) WRITE(6,'(//,A\,E12.5)') ' DET < 0  =====>>> ',DET
      CALL UPDATE(DISP)
      CALL OUTPUT

      WRITE(*,*) '*** PRESS ANY KEY ***'
      READ(*,*) 					
   
      WRITE(6,101)

      CLOSE(2);  CLOSE(4);  CLOSE(5);  CLOSE(6);

      STOP
 101  FORMAT(///12X,'*************** END OF ANALYSIS ***************'/)
 2    FORMAT(I6,2X,3E12.4)
      END

C...............................................................................

      SUBROUTINE DISPL(UR,RP)
      USE GCONTROL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RP(*),UR(*)

      WRITE(6,102)
      DO I=1,NEL+1
        WRITE(6,'(I5,4E15.6)') I,RP(I),UR(I)
      ENDDO

      WRITE(2,'(5E15.6)') (UR(I), I=1,NEL+1)
      WRITE(6,'(//,A,I6)') ' NUMBER OF DOF = ', NEL+1

  102 FORMAT(//,1X,'N O D A L    D I S P L A C E M E N T S.........',
     +  //,1X,'NODE   R-COORDINATE   DISPLACEMENT',/)
      END SUBROUTINE DISPL 

C...............................................................................

      SUBROUTINE INPUT(RG,ICO,IX)
      USE GCONTROL;USE MPROPERTY; USE MPOINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RG(*),ICO(3,*),IX(*),LJ(3)
      CHARACTER *64 TITLE

C** READ IN THE MASTER CONTROL INFORMATION

      READ(5,'(A)')  TITLE      
      READ(5,*)  NEL, IFix, Rmin, Rmax, Tr

      WRITE(6,'(/,1X,A)') TITLE
      WRITE(6,100)

C**DETERMINE NODE LOCATIONS AND ELEMENT NUMBERING, CALC. HALF BANDWIDTH

      CALL LAYOUT(RG,ICO,IX)
      CALL BANDWH(ICO,IX,LJ)
      NVA = (3*LBAND+1)*NNET

      RP(1) = Rmin;       
      DO IEL=1,NEL
	  RP(IEL+1) = RG(2*IEL+1)   ! NODES FOR DISPLACEMENT DOF
        RC(IEL) = RG(2*IEL)       ! NODES FOR PRESSURE DOF
      ENDDO

      WRITE(6,102) NEL,NNOD,NVAR,NNODEL,IFix,Rmin,Rmax,DL,TR

      WRITE(6,110) NVEL,LBAND,NVA

C**  MATERIAL PROPERY SETS

      READ(5,*) EMOD,ANV
      GMOD = EMOD/(1.D0+ANV)/2.D0
      xkap = 2.D0*GMOD
      CF = 0.D0

      If(ANV <0.4995) then
        WRITE(*,'(A,\)') ' COMPRESSIBILITY FACTOR OF FLUID = '
        READ(*,*) CF	
        XKAP = 2.D0*GMOD/(1.D0-2.D0*ANV)
      ENDIF

      WRITE(6,115); I = 1;
      WRITE(6,116) I,EMOD,GMOD,ANV,CF

      CF = 0.D0
      WRITE(6,'(//,A,E10.5)') ' COMPRESSIBILITY OF FLUID = ', CF
      Write(6,'(   A)')       ' ************************'
  
  100 FORMAT(/,5x,'***** UP-FORMULATION: AXI-SYMMETRIC *****',/)    
  102 FORMAT( /,5X,'TOTAL NO. OF ELEMENTS             NEL   =',I5,
     +        /,5X,'TOTAL NO. OF NODES                NNOD  =',I5,
     +        /,5X,'VARIABLES PER NODE                NVAR  =',I5,
     +        /,5X,'NO. OF NODES PER ELEM.            NNODEL=',I5,
     +        /,5X,'FIXITY OF LAST NODE               IFix  =',I5,
     +        /,5X,'BOREHOLE RADIUS                   Rmin  =',E12.4,
     +        /,5X,'DOMAIN RADIUS                     Rmax  =',E12.4,
     +        /,5X,'ELEMENT LENGTH                      DL  =',E12.4,
     +        /,5X,'PRESSUREMETER PRESSURE              Tr  =',E12.4,/)
  110 FORMAT(/,5X,'NO. OF VARIABLES PER ELEMENT      NVEL  =',I5,
     +       /,5X,'HALF BANDWIDTH EXCLUDING DIAG.    LBAND =',I5,
     +       /,5X,'SIZE OF GLOBAL STIFFNESS MATRIX   NVA   =',I8,/) 
 115  FORMAT(///,1X,'M A T E R I A L    P R O P E R T I E S ')
 116  FORMAT(//,'PROPERTIES FOR MATERIAL SET ',I5,
     +        /,'*************************** ',
     +        /,'ELASTIC MODULUS (EMOD)............= ',E12.5,
     +        /,'SHEAR MODULUS (GMOD)..............= ',E12.5,
     +        /,'POISSON RATIO (ANV)...............= ',E12.5,
     +        /,'COMPRESSIBILITY OF FLUID (CF).....= ',E12.5,/)

      END SUBROUTINE INPUT

C.........................................................................

      SUBROUTINE LAYOUT(RG,ICO,IX)
      USE GCONTROL;
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RG(*),IX(*),ICO(3,*)

      NNOD = 2*NEL + 1
      NNODEL = 3
      NVAR = 1
      NVEL=NVAR*NNODEL      
      NMAT=NVAR*NNOD

      ix(1:nnod) = 1;    
      IF(IFix > 0) ix(nnod) = 0;   !  maximum radius dof fixed

      RG(1) = Rmin;  DL = (Rmax-Rmin)/NEL; DR = DL/2;  

      DO I = 2, NNOD
       RG(I) = RG(I-1) + DR;
      ENDDO

      !---------------------------------------------------------
      !  input data for nodal connectivity for each element
      !  ico(i,j) where i-> element no. and j-> connected nodes
      !---------------------------------------------------------

      DO I = 1, NEL
        ico(1,i) = 2*i-1; ico(2,i) = 2*i+1; ico(3,i) = 2*i;
      ENDDO

      NMAT=NVAR*NNOD
      NNET=0
      DO I=1,NMAT
      IF(IX(I).GT.0) THEN
        NNET=NNET+1
        IX(I)=NNET
      ELSE
        IX(I)=0
      ENDIF
      ENDDO

      END SUBROUTINE LAYOUT

C.........................................................................

      SUBROUTINE MATRX(GSTIF,ICO,RG,IX)
      USE GCONTROL; USE MPROPERTY; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: GSTIF(*),ESTIF(3,3),RG(*),A,B
      INTEGER :: IEL,ICO(3,*),IX(*),LJ(3)
	
C** FORM GLOBAL STIFFNESS MATRIX

      GSTIF(1:NVA) = 0.D0
      DO IEL=1,NEL	
         A = RG(ICO(1,IEL)); B = RG(ICO(2,IEL));
         LJ(1:3)= IX(ICO(1:3,IEL))
         CALL STIFF(A,B,ESTIF)
         CALL MAPST(GSTIF,ESTIF,NVEL,LJ,LBAND)
      ENDDO

      RETURN
      END SUBROUTINE MATRX		   
		
C...............................................................................

      SUBROUTINE OUTPUT
      USE GCONTROL; USE MPROPERTY; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: P, Q, ST(4)
      INTEGER :: IEL

C** PRINT DISPLACEMENT

      CALL DISPL(UR,RP)

      rewind 4

C** STRESSES AT CENTROID

      WRITE(6,101)
      DO IEL=1,NEL
        ST(1) = SXX(IEL)- UW(IEL); ST(2) = SYY(IEL) - UW(IEL); !total stresses
        ST(3) = SXY(IEL);          ST(4) = SZZ(IEL) - UW(IEL);
        P = -(ST(1)+ST(2))/2.D0;      ! sr, etc are total stresses
        Q = SQRT((ST(1)-ST(3))**2/4.D0+ST(3)**2) 
        WRITE(6,104) IEL,RC(IEL),ST(1:4),UW(IEL)   !,P,Q ! Total stresses
        WRITE(4,'(6e15.6)') SXX(IEL),SYY(IEL),SXY(IEL),SZZ(IEL),UW(IEL)  
      ENDDO					                     

      RETURN
  101 FORMAT(//,' MATERIAL POINT STRESSES',
     +        /,' ***********************',
     +   //,1X,'MPN    REFERENCE   LOCATION   STRESS -11   STRESS -22',
     +     '   STRESS -12   STRESS -33',//)
  104 FORMAT(I5,3X,E9.4,11X,7E13.5)

      END SUBROUTINE OUTPUT

C...............................................................................

      SUBROUTINE STIFF(a,b,S)
      USE MPROPERTY  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(3,3)

      S(1:3,1:3) = 0.D0

C     THIS SUBROUTINE WILL CALCULATE THE STIFFNESS MATRIX FOR EACH ELEMENT.  

      c1 = xkap/(b-a)**2;  c2 = log(b/a);  V = (b**2-a**2)/2.D0;

      If(ANV < 0.4995) c2 = (1-anv)*c2;

      S(1,1) = c1*(b**2*c2-(b-a)**2); S(1,2) = c1*a*b*c2; S(1,3) = a;  
      S(2,1) = c1*a*b*c2; S(2,2) = c1*(a**2*c2+(b-a)**2); S(2,3) = -b; 
      S(3,1) = a;         S(3,2) = -b;                 S(3,3) = -CF*V;

      RETURN
      END SUBROUTINE STIFF

C..........................................................................

      SUBROUTINE UPDATE(U)
      USE GCONTROL; USE MPROPERTY; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: ER,ET,U(*) 
      INTEGER :: IEL

C**   UPDATE DISPLACEMENTS AND STRESSES

      IF(IFix.EQ.0) U(NNOD) = 0.D0
	UR(1) = U(1);
      DO IEL=1,NEL
	  UR(IEL+1) = U(2*IEL+1) ! DISPLACEMENT AT NODES 
	  UW(IEL) = U(2*IEL)   ! UW AT ELEMENT CENTROID
        ER = (UR(IEL+1) - UR(IEL))/DL
        ET = (UR(IEL+1) + UR(IEL))/(2*RC(IEL))

        If(ANV <0.4995) then
          SXX(IEL) = XKAP*((1.D0-ANV)*ER+ANV*ET)
          SYY(IEL) = XKAP*(ANV*ER+(1.D0-ANV)*ET)
          SZZ(IEL) = ANV*(SXX(IEL)+SYY(IEL))
          SXY(IEL) = 0.D0
        else
          SXX(IEL) = XKAP*ER
          SYY(IEL) = XKAP*ET
          SZZ(IEL) = ANV*(SXX(IEL)+SYY(IEL))
          SXY(IEL) = 0.D0
        endif
      ENDDO

      RETURN
      END SUBROUTINE UPDATE
	  
C********************************************************************
   
      SUBROUTINE GBAND(A,B,LT,IP,DET,NCN,ML,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION  B(1),IP(1),A(1)
      IFN(I,J)=1+(J-1)*LC+I-J+NUM
      NU=ML
      LCM=NU+2*ML
      LC=LCM+1
      NLC=N*LC
      NUM=NU+ML
C  GENERATE REMAINING ELEMENTS OF SYMMETRIC MATRIX
      IF(LT.NE.-1) GO TO  4
      NN=N-1
      DO 32 I=1,NN
      IFI=IFN(I,I)
      IFJ=IFI
      II=I+1
      IML=MIN0(I+ML,N)
      DO 33 J=II,IML
      IFI=IFI+1
      IFJ=IFJ+LCM
33    A(IFJ)=A(IFI)
32    CONTINUE
4     IF(IABS(LT).NE.1) GO TO 10
      IP(N)=1
      IF(ML.EQ.0) GO TO 35
C  SET ELEMENTS 1 - ML OF EACH COLUMN TO ZERO
      DO 31 I=1,N
      IFK=(I-1)*LC
      DO 31 J=1,ML
      IFK=IFK+1
31    A(IFK)=0.D0
   35 DET=0.D0
      NCN=0
      IF(ML.EQ.0) GO TO 15
C  LU DECOMPOSITION
      DO 1 K=1,N
      IFK=IFN(K,K)
      IF(K.EQ.N) GO TO 13
      KP=K+1
      KPM=MIN0(K+ML,N)
      KPN=MIN0(K+NUM,N)
      M=K
      IFM=IFK
      IFI=IFK
      DO 3 I=KP,KPM
      IFI=IFI+1
      IF(DABS(A(IFI)).LE.DABS(A(IFM))) GO TO 3
      M=I
      IFM=IFI
3     CONTINUE
      IP(K)=M
      T=A(IFM)
      IF(M.NE.K) IP(N)=-IP(N)
      A(IFM)=A(IFK)
      A(IFK)=T
      IF(T.EQ.0.0) GO TO 20
      OT=1.D0/T
      IK=IFK
      DO 5 I=KP,KPM
      IK=IK+1
5     A(IK)=-A(IK)*OT
      KJ=IFK
      MJ=IFM
      DO 6 J=KP,KPN
      KJ=KJ+LCM
      MJ=MJ+LCM
      T=A(MJ)
      A(MJ)=A(KJ)
      A(KJ)=T
      IF(T.EQ.0.D0) GO TO 6
      IK=IFK
      IJ=KJ
      DO 7 I=KP,KPM
      IK=IK+1
      IJ=IJ+1
7     A(IJ)=A(IJ)+A(IK)*T
6     CONTINUE
   13 IF(A(IFK).EQ.0) GO TO 20
1     CONTINUE
15    IFK=IFN(1,1)
      DET=A(IFK)
      DO 16 K=2,N
      IFK=IFK+LC
      DET=DET*A(IFK)
      IF(DET.EQ.0.D0) GO TO 20
      IF(DABS(DET).GT.1.D-15) GO TO 17
      DET=DET*1.D+15
      NCN=NCN-15
      GO TO 16
17    IF(DABS(DET).LT.1.D+15) GO TO 16
      DET=DET*1.D-15
      NCN=NCN+15
16    CONTINUE
      DET=DET*IP(N)
      GO TO 10
20    DET=0.D0
      WRITE(6,100) K
100   FORMAT('***MATRIX IS SINGULAR,',/,5X,
     1'ERROR OCCURRED IN ATTEMPT TO FIND',I5,'TH PIVOT')
      RETURN
10    CALL SOLVE(A,B,IP,ML,N)
      RETURN
      END SUBROUTINE GBAND

C...................................................................

      SUBROUTINE SOLVE(A,B,IP,ML,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(1),B(1),IP(1)
C  MATRIX A INTO A LOWER TRIANGULAR MATRIX L AND AN UPPER TRIANGULAR
C  MATRIX U.
C
      IFN(I,J)=1+(J-1)*LC+I-J+NUM
      NU=ML
      LCM=2*ML+NU
      LC=LCM+1
      NUM=NU+ML
      MN=N-1
C  SOLVE FOR Y
      IF(ML.EQ.0) GO TO 10
      DO 1 K=1,MN
      KP=K+1
      M=IP(K)
      T=B(M)
      B(M)=B(K)
      B(K)=T
      KPM=MIN0(K+ML,N)
      IFK=IFN(K,K)
      DO 1 I=KP,KPM
      IFK=IFK+1
1     B(I)=B(I)+A(IFK)*T
C  SOLVE FOR X
10    IFK=IFN(N,N)
      DO 2 KB=1,MN
      KM=N-KB
      K=KM+1
      B(K)=B(K)/A(IFK)
      IFK=IFK-LC
      T=-B(K)
      KMN=MAX0(1,K-ML-NU)
      KML=IFN(KMN,K)
      DO 2 I=KMN,KM
      B(I)=B(I)+A(KML)*T
2     KML=KML+1
      B(1)=B(1)/A(NUM+1)
      RETURN
      END SUBROUTINE SOLVE

C......................................................................

      SUBROUTINE BANDWH(ICO,JX,LJ)
      USE GCONTROL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ICO(3,*),JX(*),LJ(*)
      LBAND=0
      NV2=2*NVAR
      DO 3 I=1,NEL
      DO 4 J=1,NVAR
      DO 4 K=1,NNODEL
      K1=(K-1)*NVAR
      LJ(J+K1)=JX(NVAR*ICO(K,I)-NVAR+J)
    4 CONTINUE
      MAX=0
      MIN=1000000
      NV3=NVAR*NNODEL
      DO 8 J=1,NV3
      IF(LJ(J).EQ.0) GO TO 8
      IF(LJ(J)-MAX) 6,6,5
    5 MAX=LJ(J)
    6 IF(LJ(J)-MIN) 7,8,8
    7 MIN=LJ(J)
    8 CONTINUE
      NB1=MAX-MIN
      IF(NB1.GT.LBAND) LBAND=NB1
    3 CONTINUE

      RETURN
      END SUBROUTINE BANDWH

C.........................................................................

      SUBROUTINE MAPLD(B,FL,LJ,NVEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(*),FL(*),LJ(*)
      DO 12 I=1,NVEL
      LJR=LJ(I)
      IF(LJR.EQ.0) GO TO 12
      B(LJR)=B(LJR)+FL(I)
   12 CONTINUE
      RETURN
      END SUBROUTINE MAPLD

C...............................................................................

      SUBROUTINE MAPST(A,S,NVEL,LJ,LBAND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),S(NVEL,NVEL),LJ(*)
      DO 12 I=1,NVEL
       LJR=LJ(I)
       IF(LJR.LE.0) GO TO 12
        DO 11 J=1,NVEL
         LJC=LJ(J)
         IF(LJC.LE.0) GO TO 11
          L=3*LBAND*LJC+LJR-LBAND
          A(L)=A(L)+S(I,J)
   11 CONTINUE
   12 CONTINUE
      RETURN
      END SUBROUTINE MAPST

C......................................................................

      SUBROUTINE MULT1(X,Y,S,Z,M1,M2,M3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C     MULTIPLIES THE MATRICES  Y (TRANSPOSE) * X * Y .

      DIMENSION X(M1,M1),Y(M1,M2),Z(M3,M2),S(M2,M2)
      DO 1 I=1,M1
      DO 2 K=1,M2
      XX=0.D0
      DO 3 J=1,M1
    3 XX=XX+X(I,J)*Y(J,K)
    2 Z(I,K)=XX
    1 CONTINUE
      DO 4 I=1,M2
      DO 5 K=I,M2
      XX=0.D0
      DO 6 J=1,M1
    6 XX=XX+Y(J,I)*Z(J,K)
      S(I,K)=XX
    5 S(K,I)=XX
    4 CONTINUE
      RETURN
      END SUBROUTINE MULT1

C..................................................................

      SUBROUTINE GLOLOC(SXX,SYY,SZZ,SXY,STRES,INTGLO,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SXX(*),SYY(*),SZZ(*),SXY(*),STRES(4)

      IF(N.EQ.0) THEN
         SXX(INTGLO) = STRES(1)
         SYY(INTGLO) = STRES(2)
         SXY(INTGLO) = STRES(3)
         SZZ(INTGLO) = STRES(4)
      ELSE
         STRES(1) = SXX(INTGLO)
         STRES(2) = SYY(INTGLO)
         STRES(3) = SXY(INTGLO)
         STRES(4) = SZZ(INTGLO)
      ENDIF

      RETURN
      END
