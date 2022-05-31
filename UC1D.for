      MODULE GCONTROL
       INTEGER, PARAMETER :: NNG = 1000, NELT  = 1000
       INTEGER, PARAMETER :: NFD = 2000, MSTIF = 100000 
       INTEGER :: NEL,NNOD,NVAR,NNODEL,NVEL,LBAND,NMAT,NVA
       INTEGER :: NSTEP,NITER,IFix
       INTEGER :: idisp = 1, istr, NPT =100
       ! dtime = time step, TR = applied pressure
       ! Rmin = location of zero excess pressure
       ! Rmin, Rmax = locations of zero excess pressure and doimain radius
	 DOUBLE PRECISION :: Rmin,Rmax,DL,TR, Cvmax,time,dtime  
       CHARACTER (LEN = 5) :: INFILE
      END MODULE GCONTROL

	MODULE MPROPERTY
        DOUBLE PRECISION :: GMOD,EMOD,ANV,XKAP
        DOUBLE PRECISION :: Kw,Cv,BULK,Em
        DOUBLE PRECISION :: GRw = 10
	END MODULE MPROPERTY

      MODULE MPOINT
        INTEGER, PARAMETER :: IPNT=1000
        DOUBLE PRECISION, DIMENSION (IPNT)::GAMA,FBAR,UW,RC,RP
        DOUBLE PRECISION, DIMENSION (IPNT)::SXX,SYY,SXY,SZZ,UR
      END MODULE MPOINT

C******************************************************************

      PROGRAM UPC1D
      USE GCONTROL; USE MPROPERTY; USE MPOINT
				
C       ************************************************************
C       *                                                          *
C       *      2-NODED STRESS ANALYSIS FINITE ELEMENT PROGRAM      *
C       *       CONSOLIDATION WITH EFFECTIVE STRESS CONSTANT       *
C       *           and U LINEAR. This is linear elastic           *
C       *                         (May 2022)                       *
C       ************************************************************

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION (NFD) :: TLOAD,ILOAD,dUR
      DOUBLE PRECISION :: GSTIF(MSTIF),RG(NNG),CC(MSTIF)
      DOUBLE PRECISION :: DET, RATIO, A, B, D, tmp, ST(3)
      INTEGER :: I,IEL,NCN,ICO(2,NELT),LJ(2) 

      WRITE(*,'(A,\)') ' INPUT FILE NAME = '
      READ(*,'(A)') INFILE

      OPEN(2,FILE=INFILE//'.DIS',FORM='FORMATTED')  ! Initial Displacements
      OPEN(3,FILE=INFILE//'.HIS',FORM='FORMATTED')  ! History at node istr
      OPEN(4,FILE=INFILE//'.STR',FORM='FORMATTED')  ! Stresses for consolidation
      OPEN(5,FILE=INFILE//'.DAT',FORM='FORMATTED')  ! Basic FEM input
      OPEN(6,FILE=INFILE//'.OUT',FORM='FORMATTED')  ! Output at end of analysis
      OPEN(7,FILE=INFILE//'.PLT',FORM='FORMATTED')  ! Plotting output 

C** INPUT THE INFORMATION NEEDED FOR THE REQUESTED FINITE ELEMENT ANALYSIS

      CALL INPUT(RG,ICO)

C*******************************************************************

C**  STEP 1: INITIALIZES ARRAYS INCLUDING INITIAL STRESSES

      TLOAD(2:NNOD) = 0.D0;  TLOAD(1) = TR*Rmin;   ! Surface load
      ILOAD(1:NNOD) = 0.D0;  UW(1:NEL) = 0.D0;
      SXX(1:NEL) = 0.D0;     SYY(1:NEL) = 0.D0;     
      SZZ(1:NEL) = 0.D0;     dUR(1:NNOD) = 0.D0;

C**   Initial Displacements from undrained analysis

      READ(2,'(5E15.6)') (UR(I), I=1,NNOD)

C**   Initial effective stresses, plus pressures
      DO I =1, NEL
        READ(4,*) RC(I),SXX(I),SYY(I),tmp,SZZ(I),UW(I)	 
      ENDDO                    ! Tension +ve for Sij    

      WRITE(6,101)
      CALL OUTPUT
  
C*******************************************************************

C**   Determine the time step

      dtime = DL**2/(4*Cv)		 ! Critical time step

      WRITE(*,'(A\,E12.5,/)')   ' Critical time step  =  ', dtime
      WRITE(*,'(A\)') ' Provide dtime     =  ' 
      READ(*,*) dtime           
      WRITE(*,'(A\)') ' Provide nstep     =  '
      READ(*,*) nstep

      WRITE(6,'(/,A\,E12.5)')   ' Time step  =  ', dtime      

      WRITE(6,'(A\,I7,/)')    ' Number of steps   =  ', nstep
      WRITE(6,'(A\,E12.5,//)')' Maximum time      =  ', nstep*dtime

      WRITE(*,'(//)')                   

C**  STEP 2: SETUP MATRICES FOR MECHANICAL & SEEPAGE ANALYSES

      CALL MATRX(GSTIF,CC,ICO,RG)    ! Set up matrices

C*******************************************************************

C**   STEP 3:  LOOP FOR TIME STEPPING DURING CONSOLIDATION

      ILOAD(1:NNOD) = 0.D0;

C**   Displacements associated with effective stresses
C**   Calculate net load associated with displacement calculation

      CALL MatrixXa(ILOAD,UR,GSTIF,NNET,LBAND)

      DUR(1:NNOD) = dtime*(Tload(1:NNOD)-ILOAD(1:NNOD))

      DET=1.D-8  
      RATIO=1.D-8  
      CALL CHOLSQ(CC,dUR,NNOD,LBAND,2,RATIO,DET,NCN)
      UR(1:NNOD) = UR(1:NNOD) + dUR(1:NNOD) 
      CALL UPDATE                                           	  	               

      CALL OUTPUT

      UW(1:NEL) = 0.D0
      DO IEL=1,NEL                                             ! plot data
        ST(1) = SXX(IEL)- UW(IEL); ST(2) = SYY(IEL) - UW(IEL); ! total stresses
        ST(3) = SZZ(IEL) - UW(IEL);
        D = (UR(IEL)+UR(IEL+1))/2.D0;
        WRITE(7,'(7e13.5)') RC(IEL),D,ST(1:3),UW(IEL)  
      ENDDO               

      
      WRITE(*,*) '*** PRESS ANY KEY ***'
      READ(*,*) 					
   
      WRITE(6,102)

      CLOSE(2);  CLOSE(3);  CLOSE(4);  CLOSE(5);  CLOSE(6);  CLOSE(7);

      STOP

  101  FORMAT(//,' INITIAL MATERIAL POINT DISPLACEMENTS & STRESSES',
     +        /,' ************************************************')
  102 FORMAT(///12X,'*************** END OF ANALYSIS ***************'/)
 104  FORMAT(I5,6X,6E13.5)
      END

C.........................................................................
		   
      SUBROUTINE DISPL(UR,RP)
      USE GCONTROL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RP(*),UR(*)

      WRITE(6,102)
      DO I=1,NEL+1
        WRITE(6,'(I5,3E15.6)') I,RP(I),UR(I)
      ENDDO

  102 FORMAT(//,1X,'N O D A L    D I S P L A C E M E N T S.........',
     +  //,1X,'NODE   R-COORDINATE',
     +     1X,'DISPLACEMENT-U',/)
      END SUBROUTINE DISPL 

C............................................................

      SUBROUTINE CAPACITANCE(A,B,ECC)
      USE GCONTROL; USE MPROPERTY;
      IMPLICIT NONE
      DOUBLE PRECISION :: ECC(2,2), A , B

C** FORM LOCAL CAPACITANCE MATRIX
      ECC(1:2,1:2) = 0.D0
      ECC(1,1) = (B-A)*(B+2.D0*A)/6.D0;  
      ECC(2,2) = (B-A)*(A+2.D0*B)/6.D0;  

      RETURN
      END SUBROUTINE HEAT		   

C............................................................

      SUBROUTINE INPUT(RG,ICO)
      USE GCONTROL; USE MPROPERTY; USE MPOINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RG(*),ICO(2,*)
      CHARACTER *64 TITLE

C** READ IN THE MASTER CONTROL INFORMATION

      READ(5,'(A)')  TITLE      
      READ(5,*)  NEL, IFix, Rmin, Rmax, Tr, istr
      IF (IFix > 1)  IFix = 1
      IF (IFix < 0)  IFix = 0

      WRITE(6,'(/,1X,A)') TITLE
      WRITE(6,100)

C**DETERMINE NODE LOCATIONS AND ELEMENT NUMBERING, CALC. HALF BANDWIDTH

      CALL LAYOUT(RG,ICO)
      LBAND = 1;  NVA = (LBAND+1)*NNOD 

      RP(1) = Rmin;       
      DO IEL=1,NEL
        RP(IEL+1) =  RG(IEL+1)                 ! NODES FOR DISPLACEMENT DOF
        RC(IEL)   = (RP(IEL)+RP(IEL+1))/2.D0   ! NODES FOR PRESSURE DOF
      ENDDO

      WRITE(6,102) NEL,NNOD,IFix,ISTR,Rmin,Rmax,DL,TR 
      WRITE(6,110) NVEL,NNOD,LBAND,NVA

C**  MATERIAL PROPERY SETS

      READ(5,*) EMOD,ANV,Kw
      GMOD = EMOD/(1.D0+ANV)/2.D0
      XKAP = 2.D0*GMOD/(1.D0-2.D0*ANV)
      Em   = 2.D0*GMOD*(1- ANV)/(1.D0-2.D0*ANV)
      Bulk = EMOD/(3.d0*(1.D0-2.D0*ANV)); CF = 0.D0
      Cv	 = Em*Kw/GRw

      WRITE(6,115); I = 1;
      WRITE(6,116) I,EMOD,GMOD,ANV,Bulk,Kw,Cv
  
  100 FORMAT(/,5x,'***** UP-CONSOLIDATION: AXI-SYMMETRIC *****',/)    
  102 FORMAT( /,5X,'TOTAL NO. OF ELEMENTS             NEL  =',I5,
     +        /,5X,'TOTAL NO. OF NODES                NNOD =',I5,
     +        /,5X,'FIXITY OF LAST NODE               IFix =',I5,
     +        /,5X,'NODE FOR OUTPUT ON SCREEN         ISTR =',I5,
     +        /,5X,'BOREHOLE RADIUS                   Rmin =',E12.4,
     +        /,5X,'DOMAIN RADIUS                     Rmax =',E12.4,
     +        /,5X,'ELEMENT LENGTH                      DL =',E12.4,
     +        /,5X,'PRESSUREMETER PRESSURE              Tr =',E12.4,/)
  110 FORMAT(/,5X,'NO. OF VARIABLES PER ELEMENT      NVEL  =',I5,
     +       /,5X,'NET DOF                           NNET  =',I5,
     +       /,5X,'HALF BANDWIDTH EXCLUDING DIAG.    LBAND =',I5,
     +       /,5X,'SIZE OF GLOBAL STIFFNESS MATRIX   NVA   =',I8,/) 
 115  FORMAT(///,1X,'M A T E R I A L    P R O P E R T I E S ')
 116  FORMAT(//,'PROPERTIES FOR MATERIAL SET ',I5,
     +        /,'*************************** ',
     +        /,'ELASTIC MODULUS (EMOD)............= ',E12.5,
     +        /,'SHEAR MODULUS (GMOD)..............= ',E12.5,
     +        /,'POISSON RATIO (ANV)...............= ',E12.5,
     +        /,'BULK MODULUS (BULK).............. = ',E12.5,
     +        /,'PERMEABILITY (Kw)................ = ',E12.5,
     +        /,'COEFFICIENT OF CONSOLIDATION (Cv) = ',E12.5,/)

      END SUBROUTINE INPUT

C.........................................................................

      SUBROUTINE LAYOUT(RG,ICO)
      USE GCONTROL;
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RG(*),ICO(2,*)

      NNOD = NEL + 1
      NNODEL = 2
      NVAR = 1
      NVEL=NVAR*NNODEL      
      NMAT=NVAR*NNOD
      ix(1:nnod) = 1;    
      IF(IFix > 0) ix(nnod) = 0;   !  maximum radius dof fixed

      RG(1) = Rmin;  DL = (Rmax-Rmin)/NEL;  

      DO I = 2, NNOD
       RG(I) = RG(I-1) + DL;
      ENDDO

      !---------------------------------------------------------
      !  input data for nodal connectivity for each element
      !  ico(i,j) where i-> element no. and j-> connected nodes
      !---------------------------------------------------------

      DO I = 1, NEL
        ico(1,i) = i; ico(2,i) = i+1;
      ENDDO

      NMAT=NVAR*NNOD

      END SUBROUTINE LAYOUT

C.........................................................................
		   			 
      SUBROUTINE MATRX(GSTIF,CC,ICO,RG)
      USE GCONTROL; USE MPROPERTY; USE MPOINT;
      IMPLICIT NONE
      DOUBLE PRECISION :: GSTIF(*),CC(*),RG(*), A, B, DET, RATIO
      DOUBLE PRECISION :: ECC(2,2),F(NNG),ESTIF(2,2)
      INTEGER :: IEL, ICO(2,*),LJ(2),NCN

C** FORM GLOBAL STIFFNESS MATRIX
	
      GSTIF(1:NVA) = 0.D0; CC(1:NVA) = 0.D0;  F(1:NNOD) = 0.D0

C**   Element  stiffness matrix H; Lumped Capacitance Matrix ECC 

      DO IEL=1,NEL
         A = RG(IEL); B = RG(IEL+1);
         LJ(1:2)= ICO(1:2,IEL)
         CALL STIFF(A,B,ESTIF)
         CALL CAPACITANCE(A,B,ECC)
         CALL MAPST(GSTIF,CC,ECC,ESTIF,NVEL,LJ,LBAND)
      ENDDO

      CC(1:NVA) =  CC(1:NVA) + dtime*GSTIF(1:NVA)

      DET=1.D-8
      RATIO=1.D-8
      CALL CHOLSQ(CC,F,NNOD,LBAND,1,RATIO,DET,NCN)
      IF(DET.LE.0.D0) WRITE(6,'(//,A\,E12.5)') ' DET < 0  =====>>> ',DET
	
      END SUBROUTINE MATRX		   
	
C............................................................

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
        ST(3) = SZZ(IEL) - UW(IEL);
        P = -(ST(1)+ST(2))/2.D0;                   ! sr, etc are total stresses
        Q = SQRT((ST(1)-ST(3))**2/4.D0+ST(3)**2) 
        WRITE(6,104) IEL,RC(IEL),ST(1:3),UW(IEL)   !,P,Q ! Total stresses
      ENDDO					                     

      RETURN
  101 FORMAT(//,' MATERIAL POINT STRESSES',
     +        /,' ***********************',
     +   //,1X,'MPN    REFERENCE   LOCATION   STRESS -11   STRESS -22',
     +     '   STRESS -33   PRESSURE  ',//)
  104 FORMAT(I5,3X,E9.4,11X,7E13.5)

      END SUBROUTINE OUTPUT

C............................................................

      SUBROUTINE STIFF(a,b,S)
      USE MPROPERTY  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(2,2)

      S(1:2,1:2) = 0.D0

C     THIS SUBROUTINE CALCULATES STIFFNESS MATRIX FOR EACH ELEMENT.  

      c1 = xkap/(b-a)**2;  c2 = (1-anv)*log(b/a);

      S(1,1) = c1*(b**2*c2-(b-a)**2); S(1,2) = c1*a*b*c2; 
      S(2,1) = c1*a*b*c2; S(2,2) = c1*(a**2*c2+(b-a)**2);
	
      RETURN
      END SUBROUTINE STIFF

C..........................................................................
	        
      SUBROUTINE UPDATE
      USE GCONTROL; USE MPROPERTY; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: ER,ET
      INTEGER :: IEL

C**   UPDATE DISPLACEMENTS AND STRESSES

      DO IEL=1,NEL
        ER = (UR(IEL+1) - UR(IEL))/DL
        ET = (UR(IEL+1) + UR(IEL))/(2*RC(IEL))

        SXX(IEL) = XKAP*((1.D0-ANV)*ER+ANV*ET)
        SYY(IEL) = XKAP*(ANV*ER+(1.D0-ANV)*ET)
        SZZ(IEL) = ANV*(SXX(IEL)+SYY(IEL))
      ENDDO		! IEL						

      RETURN
      END SUBROUTINE UPDATE
	  
C********************************************************************   

C      UTILITY ROUTINES

C********************************************************************  

      
      SUBROUTINE CHOLSQ(A,B,NDEG,LBAND,LT,RATIO,DET,NCN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(*)
C  CHOLASKI DECOMPOSITION OF THE MATRIX A BY NOT SQUARE ROOTING.
      IF(LT.GT.1) GO TO 8
      BIG=0.D0
      SMALL=1.E+16
      DO 1 I=1,NDEG
      II=I-1
      NR=I+LBAND
      IF(NR.GT.NDEG) NR=NDEG
      DO 2 J=I,NR
      KK=(I-1)*LBAND+J
      IF(I.NE.J) GO TO 3
      SII=0.D0
      IF(II.LE.0) GO TO 5
      J1=I-LBAND
      IF(J1.LE.0) J1=1
      DO 4 IR=J1,II
      IRI=(IR-1)*LBAND+I
      IRR=(IR-1)*LBAND+IR
    4 SII=SII+A(IRI)*A(IRI)*A(IRR)
    5 A(KK)=A(KK)-SII
      IF(DABS(A(KK)).GT.BIG) BIG=DABS(A(KK))
      IF(DABS(A(KK)).LT.SMALL) SMALL=DABS(A(KK))
      GO TO 2
    3 III=(I-1)*LBAND+I
      SII=0.D0
      IF(II.LE.0) GO TO 6
      J1=J-LBAND
      IF(J1.LE.0) J1=1
      IF(J1.GT.II) GO TO 6
      DO 7 IR=J1,II
      IRI=(IR-1)*LBAND+I
      IRJ=(IR-1)*LBAND+J
      IRR=(IR-1)*LBAND+IR
    7 SII=SII+A(IRI)*A(IRJ)*A(IRR)
    6 IF(A(III).EQ.0.D0) NROW=I
      IF(A(III).EQ.0.D0) GO TO 11
      A(KK)=(A(KK)-SII)/A(III)
    2 CONTINUE
    1 CONTINUE
      RAT=SMALL/BIG
      IF(RAT.LT.RATIO) WRITE(6,14) RAT,RATIO
      IF(RAT.LT.RATIO) STOP
    8 IF(LT.EQ.0) GO TO 13
      CALL FORBAK(A,B,NDEG,LBAND)
      IF(LT.GT.1) RETURN
C  CALCULATION OF THE DETERMINANT OF MATRIX A.
   13 DET=1.D0
      NCN=0
      DO 9 I=1,NDEG
      KK=(I-1)*LBAND+I
      DET=DET*A(KK)
      IF(DABS(DET).GT.1.E-15) GO TO 10
      DET=DET*1.E+15
      NCN=NCN-15
      GO TO 9
   10 IF(DABS(DET).LT.1.E+15) GO TO 9
      DET=DET*1.E-15
      NCN=NCN+15
    9 CONTINUE
      RETURN
   11 WRITE(6,12) NROW
      DET=0.D0
      NCN=0
      IF(LT.GT.0) STOP
      RETURN
   12 FORMAT(/,5X,'MATRIX IS SINGULAR',/,5X,'ERROR CONDITION OCCURRED IN
     1ROW',I5)
   14 FORMAT(/,5X,'MATRIX IS ILL-CONDITIONED,PROGRAM HALTED.',/,
     15X,'COMPUTED RATIO RAT =',E12.4,2X,'IS LESS THAN',/,
     25X,'THE SPECIFIED RATIO =',E12.4)
      END

C.........................................................................

      SUBROUTINE FORBAK(A,B,NDEG,LBAND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(*)
C  FORWARD SUBSTITUTION FOR THE INTERMEDIATE SOLUTION VECTOR
C  AND BACKWARD SUBSTITUTION FOR THE SOLUTION VECTOR.
      DO 1 I=1,NDEG
      IRII=(I-1)*LBAND+I
      J2=I-1
      IF(J2.GT.NDEG) J2=NDEG
      J1=I-LBAND
      IF(J1.LE.0) J1=1
      SII=0.D0
      IF(J1.GT.J2) GO TO 1
      DO 2 IR=J1,J2
      IRI=(IR-1)*LBAND+I
      IRR=(IR-1)*LBAND+IR
    2 SII=SII+A(IRI)*A(IRR)*B(IR)
    1 B(I)=(B(I)-SII)/A(IRII)
      DO 3 I=1,NDEG
      II=NDEG-I+1
      IJ=II+1
      NR=II+LBAND
      IF(NR.GT.NDEG) NR=NDEG
      SII=0.D0
      IF(IJ.GT.NDEG) GO TO 4
      DO 5 IR=IJ,NR
      IRI=(II-1)*LBAND+IR
    5 SII=SII+A(IRI)*B(IR)
    4 B(II)=B(II)-SII
    3 CONTINUE
      RETURN
      END
	  
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
      END

C.........................................................................
                     
      SUBROUTINE MAPST(A,C,ECC,S,NVEL,LJ,LBAND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),C(*),S(NVEL,NVEL),ECC(NVEL,NVEL),LJ(*)
      DO 12 I=1,NVEL
      LJR=LJ(I)
      IF(LJR.EQ.0) GO TO 12
      DO 11 J=I,NVEL
      LJC=LJ(J)
      IF(LJC.EQ.0)  GO TO 11
      IF(LJR-LJC) 9,10,10
   10 K=(LJC-1)*LBAND+LJR
      GO TO 13
    9 K=(LJR-1)*LBAND+LJC
   13 A(K)=A(K)+S(I,J)
      C(K)=C(K)+ECC(I,J)

   11 CONTINUE
   12 CONTINUE
      RETURN
      END

C..................................................................

      SUBROUTINE MatrixXa(B,X,S,NNET,LBAND) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(1), B(1), X(1)
      B(1:NNET) = 0.d0

C **  NOT USED FOR GLOBAL MULTIPLICATION IN THIS CODE
C **  CALL MatrixXa(BB,TDISP,GSTIF,NNET,LBAND) 
   
      DO I = 1, NNET	        !  B = S*X
        LC1 = I - LBAND
        IF(LC1.LE.0) LC1 = 1
        LC2 = I + LBAND
        IF(LC2.GT.NNET) LC2 = NNET
        B(I) = 0.D0
        DO J =  LC1, LC2
          IF((I-J).LT.0) THEN
            IJ = LBAND*(I-1) + J
          ELSE
            IJ = LBAND*(J-1) + I
          ENDIF
        ENDDO
        B(I) = B(I) + S(IJ)*X(J)
      ENDDO

      RETURN
      END
