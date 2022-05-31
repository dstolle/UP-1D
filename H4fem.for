      MODULE GCONTROL
       INTEGER, PARAMETER :: NNG = 4000, NELT  = 4000
       INTEGER, PARAMETER :: NFD = 9000, MSTIF = 3000000 
       INTEGER :: NEL,NNOD,NVAR,NNODEL,NVEL,NNET,NITER,LBAND,NVA
       INTEGER :: NCON,IELPLOT,IREAD,ITYPE,NLN
C       CHARACTER (LEN = 64) :: INFILE
       CHARACTER (LEN = 5) :: INFILE
      END MODULE GCONTROL

	MODULE MPROPERTY
      DOUBLE PRECISION, DIMENSION (20)::GMOD,EMOD,ANV,GRR,CM
	  DOUBLE PRECISION :: SIG0, beta_B
	END MODULE MPROPERTY

      MODULE MPOINT
        INTEGER, PARAMETER :: IPNT=16000
        DOUBLE PRECISION, DIMENSION (IPNT)::GAMA,FBAR,PRESSURE
        DOUBLE PRECISION, DIMENSION (IPNT)::SXX,SYY,SXY,SZZ,XMC,YMC
      END MODULE MPOINT

	MODULE TRACTIONS
        DOUBLE PRECISION, DIMENSION (2,100)::TNF,TSF
        INTEGER :: ICOB(2,100),NELB
	END MODULE TRACTIONS

      MODULE joint
        Double Precision, dimension (400)  ::KNN,KNS,TNN,TNS
        Integer :: ICOJ(2,400),NELJ
      END MODULE joint

	MODULE NUM_INT
        DOUBLE PRECISION, DIMENSION (4)::AN,ANS,ANT
        DOUBLE PRECISION, DIMENSION (4)::WF,XS,XT
        DOUBLE PRECISION :: AREA(8000)
        INTEGER :: NINTP = 4
	END MODULE NUM_INT

C******************************************************************

      PROGRAM HM4fem
      USE GCONTROL; USE MPROPERTY; USE NUM_INT
      USE TRACTIONS; USE MPOINT; USE JOINT
				
C       ************************************************************
C       *                                                          *
C       *      4-NODED STRESS ANALYSIS FINITE ELEMENT PROGRAM      *
C       *              HYBRID ELEMENT WITH SIJ LINEAR              *
C       *                    This is linear elastic                *
C       *                       (December 2021)                    *
C       ************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION (NFD) ::  TLOAD,GLOAD,DISP,TDISP
      DOUBLE PRECISION :: XG(NNG),YG(NNG)
      DOUBLE PRECISION :: GSTIF(MSTIF)
      DOUBLE PRECISION :: DET, RATIO
      INTEGER :: NCN, ismooth
      INTEGER :: ico(5,NELT),IX(NFD)

      WRITE(*,'(A,\)') ' INPUT FILE NAME = '
      READ(*,'(A)') INFILE

      OPEN(4,FILE=INFILE//'.HIS',FORM='FORMATTED')
      OPEN(5,FILE=INFILE//'.DAT',FORM='FORMATTED')
      OPEN(6,FILE=INFILE//'.OUT',FORM='FORMATTED')
      OPEN(7,FILE=INFILE//'.NND',FORM='FORMATTED')
      OPEN(8,FILE=INFILE//'.IEL',FORM='FORMATTED')
      OPEN(9,FILE = 'smooth.iel',FORM='FORMATTED')
 
	GOTO 500
      OPEN(3,FILE='INFEM.TXT',FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  PLT
C        OPEN(2,FILE=INFILE,FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  DAT
        OPEN(5,FILE=INFILE,FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  OUT
        OPEN(6,FILE=INFILE,FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  HIS
        OPEN(4,FILE=INFILE,FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  STN
C        OPEN(4,FILE=INFILE,FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  NND
        OPEN(7,FILE=INFILE,FORM='FORMATTED')
	  READ(3,'(A)') INFILE						!  IEL
        OPEN(8,FILE=INFILE,FORM='FORMATTED')
	CLOSE(3)
  500 CONTINUE

C** INPUT THE INFORMATION NEEDED FOR THE REQUESTED FINITE ELEMENT ANALYSIS

      CALL INPUT(XG,YG,ICO,IX)
      CALL GAUSS
	
      IF(NVA > MSTIF .OR. NNOD > NNG .OR. NEL > NELT) THEN
        WRITE(*,*) '*** YOU HAVE EXCEEDED SIZE LIMITS, PRESS A KEY ***'
        READ(*,*) 					
        STOP
      ENDIF  

C** INITIALIZES MOST ARRAYS INCLUDING INITIAL STRESSES

      TDISP(1:NNET) = 0.D0;  DISP(1:NNET)  = 0.D0

      CALL LOAD(GLOAD,TLOAD,XG,YG,ICO,IX)
      
C**   LOOP FOR NONLINEAR STEPPING

      WRITE(*,'(//)')


      TDISP(1:nnet) = TLOAD(1:nnet) + GLOAD(1:nnet)


      DET=1.D-8
      RATIO=1.D-8
      CALL MATRX(GSTIF,ICO,XG,YG,IX)
      CALL CHOLSQ(GSTIF,TDISP,NNET,LBAND,1,RATIO,DET,NCN)
      IF(DET.LE.0.D0) WRITE(6,'(//,A\,E12.5)') ' DET < 0  =====>>> ',DET
      CALL UPDATE(XG,YG,ICO,IX,TDISP)

      WRITE(*,'(A,\)') ' Smoothen output (ismooth = 1) = '
      READ(*,*) ismooth

      IF(ismooth ==1) then  
        CALL SMOOTH(ICO,XG,YG,TDISP,IX)
      ELSE
        CALL PLOT(XG,YG,ICO,IX,TDISP)
      ENDIF
      CALL OUTPUT(XG,YG,ICO,TDISP,IX)

      WRITE(*,*) '*** PRESS ANY KEY ***'
      READ(*,*) 					
   
      WRITE(6,101)

      CLOSE(4);  
      CLOSE(5);  CLOSE(6);  CLOSE(7);  CLOSE(8); CLOSE(9)

      STOP
 101  FORMAT(///12X,'*************** END OF ANALYSIS ***************'/)
 2    FORMAT(I6,2X,3E12.4)
      END

C.............................................................

      SUBROUTINE BMATRX(X,Y,B,DET)
      USE GCONTROL; USE NUM_INT
      IMPLICIT NONE
      DOUBLE PRECISION :: X(4),Y(4),B(3,*)
      DOUBLE PRECISION :: J11,J12,J21,J22, DET, AI(2,2)

      B(1:3,1:8) = 0.D0
      J11 = DOT_PRODUCT(ANS,X); J12 = DOT_PRODUCT(ANS,Y)
      J21 = DOT_PRODUCT(ANT,X); J22 = DOT_PRODUCT(ANT,Y)

      DET = J11*J22 - J12*J21 
      AI(1,1) =  J22/DET;  AI(1,2) = -J12/DET                   
      AI(2,1) = -J21/DET;  AI(2,2) =  J11/DET                   
      
      B(1,1:7:2)=AI(1,1)*ANS(1:4)+AI(1,2)*ANT(1:4)       
      B(3,1:7:2)=AI(2,1)*ANS(1:4)+AI(2,2)*ANT(1:4)       

      B(2,2:8:2)=B(3,1:7:2)  
      B(3,2:8:2)=B(1,1:7:2)

      RETURN
      END SUBROUTINE BMATRX 

C.......................................................................

      SUBROUTINE bmatrxJ(X,Y,B,XL)
      USE gcontrol;
      Implicit None
      DOUBLE PRECISION :: F,X(*),Y(*),B(2,4),XL,S,C

      F =  -0.5d0
      B(1:2,1:4) = 0.d0

      XL = sqrt((x(2)-x(1))**2+(y(2)-y(1))**2)
      C = (x(2)-x(1))/XL
      S = (y(2)-y(1))/XL
      
      B(1,1) = c*F;  B(1,2) = s*F; B(1,3) = c*F;  B(1,4) = -s*F;
      B(2,1) = -s*F; B(2,2) = c*F; B(2,3) = -s*F; B(2,4) =  c*F;
	 
      END	SUBROUTINE bmatrxJ


C.........................................................

      SUBROUTINE BONDRY(FL,X,Y,TN1,TS1)
      USE GCONTROL
      IMPLICIT NONE
      DOUBLE PRECISION :: FL(*),X(*),Y(*),EN(2),ENS(2),S(2),W(2)
      DOUBLE PRECISION :: TN,TS,TN1,TS1,TN2,TS2,TX,TY,DXS,DYS 
      INTEGER :: I, NODE = 2

      S = (/-0.577350269190, 0.577350269190/)
      W = (/ 1.D0, 1.D0/)

      FL(1:8) = 0.D0
      TN2 = TN1;  TS2 = TS1

      DO I = 1, NODE
        CALL ESHAPE(EN,ENS,S(I))
        DXS=ENS(1)*X(1)+ENS(2)*X(2)
        DYS=ENS(1)*Y(1)+ENS(2)*Y(2)
        TN = EN(1)*TN1+EN(2)*TN2
        TS = EN(1)*TS1+EN(2)*TS2
        TX=DYS*TN+DXS*TS
        TY=-DXS*TN+DYS*TS

        FL(1:3:2) = FL(1:3:2) +  EN(1:2)*TX*W(I)
        FL(2:4:2) = FL(2:4:2) +  EN(1:2)*TY*W(I)

      ENDDO

      RETURN
      END SUBROUTINE BONDRY

C...............................................................................

      SUBROUTINE CENTROID(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XB,YB,AREA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      X2 = X2-X1; X3 = X3-X1; X4 = X4 - X1;
      Y2 = Y2-Y1; Y3 = Y3-Y1; Y4 = Y4 - Y1;
      T0 = X2*X3; T2 = X4*X3; T4 = X3**2; T6 = X4**2; T8 = X2**2
      T3 = Y4**2; T5 = Y2**2; T7 = Y3**2; T9 = Y3*X3
      A6 = 3.D0*(X2*Y3-X3*Y2+X3*Y4-X4*Y3)
      AREA = (-X3*Y2+X3*Y4+X2*Y3-X4*Y3)/2.D0
      XB = (-T0*Y2+T2*Y4-T2*Y3+T0*Y3+T4*Y4-T6*Y3+T8*Y3-T4*Y2)/A6+X1
      YB = (Y2*X2*Y3+X3*T3-X3*T5-X4*T7-T9*Y2-Y4*X4*Y3+T9*Y4+X2*T7)/A6+Y1
      RETURN
      END SUBROUTINE CENTROID

C........................................................

      SUBROUTINE CLOAD(B,IX,NLN,NVAR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION B(1),U(NVAR),IX(1)
      IF(NLN.EQ.0) RETURN
      DO 1 I=1,NLN
      READ(5,*) INODE,(U(J),J=1,NVAR)
      NN=NVAR*(INODE-1)
      DO 2 J=1,NVAR
      N1=IX(NN+J)
      IF(N1.EQ.0) GO TO 2
      B(N1)=B(N1)+U(J)
    2 CONTINUE
    1 CONTINUE
      RETURN
      END SUBROUTINE CLOAD

C...............................................................................

      SUBROUTINE CMATRX(gmod,anu,C)
      DOUBLE PRECISION :: gmod,anu,C(3,3)
      C(1:3,1:3) = 0.D0
      IF(anu.LE.0.4995) then
      C(1,1) = (1.d0-anu)/(2.D0*gmod); C(1,2) =       -anu/(2.D0*gmod);
      C(2,1) =       -anu/(2.D0*gmod); C(2,2) = (1.d0-anu)/(2.d0*gmod);
      ELSE
        C(1,1) = 1.D0/(2.d0*gmod)
        C(2,2) = 1.D0/(2.d0*gmod)
      ENDIF
      C(3,3) = 1.d0/gmod
      RETURN
      END SUBROUTINE CMATRX

C...............................................................................

      SUBROUTINE DISPL(DIS,IX,XG,YG)
      USE GCONTROL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DIS(*),IX(*),D(4),XG(*),YG(*)

      WRITE(6,102)
      DO 1 I=1,NNOD
      I2=NVAR*I
      I1=I2-NVAR+1
      II=1
      DO 2 J=I1,I2
      D(II)=0.D0
      IF(IX(J).NE.0) D(II)=DIS(IX(J))
      II=II+1
 2    CONTINUE

      WRITE(6,103) I,XG(I),YG(I),(D(K),K=1,NVAR)

 1    CONTINUE
      RETURN

  102 FORMAT(//,1X,'N O D A L    D I S P L A C E M E N T S.........',
     +       //,1X,'NODE   X-COORDINATE   Y-COORDINATE',
     +          1X,'DISPLACEMENT-U  DISPLACEMENT-V',/)
  103 FORMAT(I5,4E15.6)

      END

C........................................................................

      SUBROUTINE ESHAPE(EN,ENS,X)
      DOUBLE PRECISION :: X, EN(*),ENS(*)

      EN(1)=(1.D0-X)/2.D0
      EN(2)=(1.D0+X)/2.D0
      ENS(1)=-0.5D0
      ENS(2)= 0.5D0

      RETURN
      END
				 				  
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
	  
C...............................................................................

      SUBROUTINE INPUT(XG,YG,ICO,IX)
      USE GCONTROL;USE MPROPERTY;USE NUM_INT;USE TRACTIONS;
      USE MPOINT; USE JOINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XG(*),YG(*),ico(5,*),IX(*),LJ(10)
      CHARACTER *64 TITLE

      WRITE(6,100)

C** READ IN THE MASTER CONTROL INFORMATION

      READ(5,'(A)')     TITLE
      WRITE(6,'(1X,A)') TITLE
      
      READ(5,*)    NEL,NNOD,NELB,NLN,NCON,NELJ,MATR	
      WRITE(6,101) MATR,NELB,NLN,NCON,NELJ
	
C**READ IN NODE LOCATIONS AND ELEMENT NUMBERING, CALC. HALF BANDWIDTH

      CALL LAYOUT(XG,YG,ICO,IX)
      CALL BANDWH(ICO,IX,LJ)
      NVA = (LBAND+1)*NNET 

      WRITE(6,110) NNET,NVA

C**  MATERIAL PROPERY SETS

      WRITE(6,115)
      DO I=1,MATR
        READ(5,*) EMOD(I),ANV(I),GRR(I)
        GMOD(I) = EMOD(I)/(1.D0+ANV(I))/2.D0
        WRITE(6,116) I,EMOD(I),ANV(I),GRR(I)
      ENDDO
      
      IF(NELB.GT.0) THEN
        WRITE(6,118)  NELB
        DO I=1,NELB
          READ(5,*) (ICOB(J,I),J=1,2),TNF(1,I),TSF(1,I)
          WRITE(6,1)I,(ICOB(J,I),J=1,2),TNF(1,I),TSF(1,I)
        ENDDO      
      ENDIF

      IF(NELJ.GT.0) THEN
        WRITE(6,119)
        DO I=1,NELJ
          READ(5,*)   (ICOJ(J,I),J=1,2),KNN(I),KNS(I)
          WRITE(6,1)I,(ICOJ(J,I),J=1,2),KNN(I),KNS(I)
        enddo      
      ENDIF

  100 FORMAT(///,'**********  F I N I T E   E L E M E N T  '
     +              ,'A N A L Y S I S   ***********',/)
  101 FORMAT(/,1X,'MIXED HYBRID STRESS ANALYSIS ',
     +       /,1X,'****************************',
     +       /,1X,'NUMBER OF DIFFERENT MATERIALS (MATR).......= ',I7,
     +       /,1X,'NUMBER OF LOADED ELEMENTS (NELB)...........= ',I7,
     +       /,1X,'Number of loaded nodes (NLN) ...........   = ',I7,
     +       /,1X,'Number of constrained dof (ncon)...........= ',I7,
     +       /,1X,'Number of interface elements (nelj)........= ',I7)

 110  FORMAT(//,'NET DEGREES OF FREEDOM(NNET).............=',I5,
     +        /,'SIZE OF GLOBAL STIFFNESS MATRIX (NVA)....=',I8)
 115  FORMAT(///,1X,'M A T E R I A L    P R O P E R T I E S ')
 116  FORMAT(//,'PROPERTIES FOR MATERIAL SET ',I5,
     +        /,'*************************** ',
     +        /,'ELASTIC MODULUS (EMOD)............= ',E12.5,
     +        /,'POISSON RATIO (ANV)...............= ',E12.5,
     +        /,'UNIT WEIGHT (GRR)................ = ',E12.5,/)
 118  FORMAT(//,' BOUNDARY LOAD NODE NUMBERS, NELB = ',I5,
     +        /,' ************************** ',/,
     +       ' ELEMENT      ICOB',5X,'  TN(1)        TS(1)')
  119  FORMAT(//,'JOINT NODE NUMBERS AND STIFFNESSES',/,
     +          '***********************************',/,
     +       ' ELEMENT      ICOJ',10X,'   KNN          KNS')
   1  FORMAT(I5,5X,2I5,1X,4E13.5,I5)
      END SUBROUTINE INPUT

C.........................................................................

      SUBROUTINE LAYOUT(XG,YG,ICO,IX)
      USE GCONTROL; USE NUM_INT; USE MPOINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XG(*),YG(*),IX(*),ico(5,*),X(4),Y(4)

      NNODEL = 4
      NVAR = 2
      NNN=NNODEL+1
      NVEL=NVAR*NNODEL

      WRITE(6,41) NEL,NNOD,NVAR,NNODEL

C     READ IN DATA
      DO I=1,NNOD
        I2=NVAR*I
        I1=I2-1
        READ(5,*) XG(I),YG(I),(IX(J),J=I1,I2)  !,IPP(I)
      ENDDO

C     WRITE INPUT

      WRITE(6,100)
      WRITE(6,42)
      DO I=1,NNOD
        I2=NVAR*I
        I1=I2-1
        WRITE(6,44) I,XG(I),YG(I),(IX(J),J=I1,I2)
      ENDDO
      
      WRITE(6,47)
      DO IEL = 1, NEL
        READ(5,*) (ICO(J,IEL),J=1,NNN)
	  X(1:4) = XG(ICO(1:4,IEL));  Y(1:4) = YG(ICO(1:4,IEL))
        CALL CENTROID(X(1),X(2),X(3),X(4),Y(1),Y(2),Y(3),Y(4),
     +                XMC(IEL),YMC(IEL),AREA(IEL))
        WRITE(6,48) IEL,(ICO(J,IEL),J=1,NNN),AREA(IEL),XMC(IEL),YMC(IEL)
      ENDDO

      NMAT=NVAR*NNOD
      NNET=0
      DO 12 I=1,NMAT
      IF(IX(I).GT.0) THEN
        NNET=NNET+1
        IX(I)=NNET
      ELSE
        IX(I)=0
      ENDIF
   12 CONTINUE

   41 FORMAT(//,'TOTAL NO. OF ELEMENTS..........',I5,/,
     +          'NUMBER OF NODES................',I5,/,
     +          'VARIABLES PER NODE.............',I5,/,
     +          'NODES PER ELEMENT..............',I5,/)
   42 FORMAT(/,3X,' NODE',6X,' X-CORD',6X,'Y-CORD',7X,' U   V')	
   44 FORMAT(1X,I5,5X,F10.3,2X,F10.3,5X,2I4)	    
   47 FORMAT(//,5X,'ELEMENT',9X,'NODE NUMBERS',5X,'IS   AREA',//)
   48 FORMAT(5X,I5,6X,4I5,I6,4E12.4)
 100  FORMAT(//,'NODAL COORDINATES AND BOUNDARY CONDITIONS',
     +        /,'*****************************************')
      END SUBROUTINE LAYOUT

C..................................................

      SUBROUTINE Lmat_B(IEL,L_B)
      USE MPROPERTY; USE NUM_INT; USE MPOINT      
      IMPLICIT NONE
      DOUBLE PRECISION :: L_B(8)
      INTEGER :: IEL
	  
      L_B(1:7:2) = 0.D0;   L_B(2:8:2) = -Area(IEL)

      RETURN
      END SUBROUTINE Lmat_B	  
	  
C...............................................................................

      SUBROUTINE LOAD(GLOAD,TLOAD,XG,YG,ICO,IX)
      USE GCONTROL; USE MPROPERTY; USE TRACTIONS
      USE MPOINT; USE NUM_INT
      IMPLICIT NONE
      DOUBLE PRECISION :: TLOAD(*),XG(*),YG(*),FL(10)
      DOUBLE PRECISION :: GLOAD(*),X(5),Y(5),DET,B(3,8)
      INTEGER :: I,IEL, IS, ico(5,*),IX(*),LJ(10), TEST
	  
      TEST = 0

      GLOAD(1:NNET) = 0.D0;  TLOAD(1:NNET) = 0.D0
      SXX(1:NEL) = 0.D0;  SYY(1:NEL) = 0.D0;     
      SXY(1:NEL) = 0.D0;  SZZ(1:NEL) = 0.D0;
      PRESSURE(1:NEL) = 0.D0;	     

C** FORM LOADS DUE TO TRACTIONS BOUNDARY TRACTIONS

      TLOAD(1:NNET) = 0.D0
      IF(NELB.GT.0) THEN
        DO I=1,NELB
          CALL LOCALB(XG,YG,X,Y,I,LJ,IX)
          CALL BONDRY(FL,X,Y,TNF(1,I),TSF(1,I))
          CALL MAPLD(TLOAD,FL,LJ,4)
        ENDDO
      ENDIF

C** FORM LOADS DUE TO POINT LOADS

      CALL CLOAD (TLOAD,IX,NLN,NVAR)
	
      RETURN
       
      IF(TEST > 0) THEN   !  This part is skipped at this time
	  
      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC(IEL),YMC(IEL))
        FL(1:NVEL) = 0.D0;
        
        DO I=1,NINTP
          CALL SHAPE(XS(I),XT(I))
          CALL BMATRX(X,Y,B,DET)          
          DET=GRR(IS)*DET
          FL(2:8:2) = FL(2:8:2)-DET*WF(I)*AN(1:4)
        ENDDO

        CALL MAPLD(GLOAD,FL,LJ,NVEL)
      ENDDO
      ENDIF

!     THE GRAVITY ROUTINES MUST BE REDONE

c      DO IEL=1,NEL
c        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC(IEL),YMC(IEL))
c        FL(1:NVEL) = 0.D0;
c        FL(2:8:2) = - GRR(IS)*AREA(IEL)/4.D0	    
c        CALL MAPLD(GLOAD,FL,LJ,NVEL)
c      ENDDO

C     WRITE THE GRAVITY LOAD VECTOR
C
C       WRITE(6,110)
C       WRITE(6,95) (GLOAD(I),I=1,NNET)


C     WRITE THE GLOBAL LOAD VECTOR

C      WRITE(6,111)
C      WRITE(6,96) (TLOAD(I),I=1,NNET)
	
   95 FORMAT(2X,4D15.5)
   96 FORMAT(2X,4D15.5)
  110 FORMAT(//,5X,'GRAVITY LOAD VECTOR',/)
  111 FORMAT(//,5X,'GLOBAL LOAD VECTOR',/)
      END SUBROUTINE LOAD

C.........................................................................

      SUBROUTINE MATRX(GSTIF,ICO,XG,YG,IX)
      USE GCONTROL; USE MPROPERTY; USE MPOINT; USE JOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: GSTIF(*),ESTIF(8,8),XG(*),YG(*)
      DOUBLE PRECISION :: X(5),Y(5),ESTIFJ(4,4),ELOAD(10)
      INTEGER :: IEL, IS, ico(5,*),IX(*),LJ(10)

C** FORM GLOBAL STIFFNESS MATRIX
	
      GSTIF(1:NVA) = 0.D0
      DO IEL=1,NEL
         CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC(IEL),YMC(IEL))
         CALL STIFF_GRAVITY(IS,X,Y,ESTIF,ELOAD)
         CALL MAPST(GSTIF,ESTIF,NVEL,LJ,LBAND)
C         CALL MAPLD(GLOAD,ELOAD,LJ,NVEL)
      ENDDO

      IF(NELJ > 0) then
        DO IEL=1,NELJ
          CALL LOCALJ(XG,YG,X,Y,IEL,LJ,IX)         
          CALL STIFFJ(IEL,X,Y,ESTIFJ)
          CALL MAPST(GSTIF,ESTIFJ,4,LJ,LBAND)
        ENDDO
      ENDIF
	
      RETURN
      END SUBROUTINE MATRX		   
	
C...............................................................................

      SUBROUTINE OUTPUT(XG,YG,ICO,TDISP,IX)
      USE GCONTROL; USE MPROPERTY; USE NUM_INT; USE MPOINT; USE JOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: X(6),Y(6),XG(*),YG(*),XXX,YYY
      DOUBLE PRECISION :: TDISP(*),ST(4),PT
      INTEGER :: IEL,IS
      INTEGER :: IX(*),ico(5,*), LJ(12)

C** PRINT DISPLACEMENT

      CALL DISPL(TDISP,IX,XG,YG)

      rewind 4

      WRITE(6,101)
      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC(IEL),YMC(IEL))
        CALL GLOLOC(SXX,SYY,SZZ,SXY,ST,IEL,1)
        WRITE(6,104) IEL,XMC(IEL),YMC(IEL),ST(1:4)
        PT = -(ST(1)+ST(2)+ST(4))/3.d0
        WRITE(4,'(5e12.4)') ST(1:4),PT
      ENDDO

C** WRITE STRESSES IN JOINT ELEMENTS

      IF(NELJ > 0) then
      WRITE(6,102)
      DO IEL=1,NELJ
        CALL LOCALJ(XG,YG,X,Y,IEL,LJ,IX)         

        xxx = (x(1)+x(2))/2   
        yyy = (y(1)+y(2))/2
 	   
        WRITE(6,104) iel,XXX,YYY,TNS(IEL),TNN(IEL)
      enddo
      ENDIF


  101 FORMAT(//,' MATERIAL POINT STRESSES',
     +        /,' ***********************',
     +   //,1X,'MPN    REFERENCE   LOCATION   STRESS -11   STRESS -22',
     +     '   STRESS -12   STRESS -33',//)
  102 FORMAT(//,' Interface Stresses',
     +        /,' ******************',
     +   //,1X,'IEL    REFERENCE   LOCATION',10X,'TNS          TNN',//)
  104 FORMAT(I5,3X,E9.3,2X,E9.3,5E13.4)
      END SUBROUTINE OUTPUT


C..........................................................................

      SUBROUTINE PLOT(XG,YG,ICO,IX,DISP)
      USE GCONTROL; USE MPROPERTY; USE NUM_INT; USE MPOINT; USE JOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: X(8),Y(8),XG(*),YG(*),V(10),V1(8),XP(4),YP(4)
      DOUBLE PRECISION :: DISP(*),DS(4),Q,B(7),C(7,8)
      DOUBLE PRECISION :: P,XL,YL,DET
      DOUBLE PRECISION :: H(7,7),Z(3,7),CMAT(3,3),Z1(3,7),R1(7,8)
      DOUBLE PRECISION :: BMAT(3,8),R(7,8),XI(3),W(3),H1(7,7)
      INTEGER :: I, IEL,J,K, L, M, IS,ico(5,1)
      INTEGER :: IX(*),LJ(18)
      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/

C**   UPDATE STRESSES

      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC(IEL),YMC(IEL))
        XP(1:4) = X(1:4) + XMC(IEL); YP(1:4) = Y(1:4) + YMC(IEL); 
        V(1:10) = 0.D0
        DO J=1,10
          IF(LJ(J) > 0) V(J)=DISP(LJ(J))
        ENDDO
        V1(1:8) = V(1:8)

        H(1:7,1:7)   = 0.D0;  H1(1:7,1:7)  = 0.D0
        R(1:7,1:8)   = 0.D0;  R1(1:7,1:8)  = 0.D0


        DO I=1,3
        DO J=1,3
          CALL SHAPE(XI(I),XI(J))
          CALL SHAPEST(X,Y,XL,YL,Z)
          CALL BMATRX(X,Y,BMAT,DET)
          CALL CMATRX(GMOD(IS),ANV(IS),CMAT)
          CALL MULT1(CMAT,Z,H1,Z1,3,7,3)

          DO 3 K=1,7
            DO 4 L=1,7
    4        H(K,L)=H(K,L)+(W(I)*W(J)*H1(K,L)*DET)
    3     CONTINUE

          DO K=1,7
            DO L = 1,8
              R1(K,L)  = DOT_PRODUCT(Z(:,K),BMAT(:,L))
            ENDDO
          ENDDO

          DO K=1,7
            DO L=1,8
              R(K,L)=R(K,L)+(W(I)*W(J)*R1(K,L)*DET)
            ENDDO
          ENDDO

        ENDDO
      ENDDO

      CALL INVERT (H,7)

        DO K=1,7
          DO L = 1,8
            C(K,L)  = DOT_PRODUCT(H(K,:),R(:,L))
          ENDDO
        ENDDO

        DO K=1,7
            B(K)  = DOT_PRODUCT(C(K,:),V1(:))
        ENDDO

        P = V(9); 

        DO M = 1,4
        DS(1) = B(1)+B(2)*X(M)+B(3)*Y(M)     
        DS(2) = B(4)+B(5)*X(M)+B(6)*Y(M)     
        DS(3) = B(7)-B(6)*X(M)-B(2)*Y(M)        
        DS(4) = ANV(IS)*(DS(1)+DS(2))
        Q = SQRT((DS(1)-DS(2))**2/4.D0+DS(3)**2)	
	  	   
        Write(7,'(7E12.4)') XP(M),YP(M),DS(1),DS(2),DS(3),P,Q

	  ENDDO  ! M LOOP
      ENDDO  ! IEL LOOP
	  
      END SUBROUTINE PLOT

C...................................................................      

      SUBROUTINE SHAPE(S,T)
      USE NUM_INT
      DOUBLE PRECISION :: S, T

      AN(1) = (1.D0-S)*(1.D0-T)/4.D0
      AN(2) = (1.D0+S)*(1.D0-T)/4.D0
      AN(3) = (1.D0+S)*(1.D0+T)/4.D0
      AN(4) = (1.D0-S)*(1.D0+T)/4.D0

      ANS(1) = -(1.D0-T)/4.D0
      ANS(2) =  (1.D0-T)/4.D0
      ANS(3) =  (1.D0+T)/4.D0
      ANS(4) = -(1.D0+T)/4.D0
      ANT(1) = -(1.D0-S)/4.D0
      ANT(2) = -(1.D0+S)/4.D0
      ANT(3) =  (1.D0+S)/4.D0
      ANT(4) =  (1.D0-S)/4.D0
      
      RETURN
      END SUBROUTINE SHAPE

C...................................................................      

      SUBROUTINE SHAPEST(X,Y,XL,YL,N)
      USE NUM_INT
      DOUBLE PRECISION :: XL,YL,X(*),Y(*),N(3,7)

!     Provides P in Matlab code; i.e., Sig = P*beta
      XL = AN(1)*X(1)+AN(2)*X(2)+AN(3)*X(3)+AN(4)*X(4)
      YL = AN(1)*Y(1)+AN(2)*Y(2)+AN(3)*Y(3)+AN(4)*Y(4)

      N(1:3,1:7) = 0.D0
      N(1,1)=1.D0; N(1,2) = XL;   N(1,3) = YL 
      N(2,4)=1.D0; N(2,5) = XL;   N(2,6) = YL 
      N(3,2)=-YL ; N(3,6) = -XL;  N(3,7) = 1.D0 
      
      RETURN
      END SUBROUTINE SHAPEST

C...............................................................
     
      SUBROUTINE SMOOTH(ICO,XG,YG,DISP,IX)
      USE GCONTROL; USE MPROPERTY; USE MPOINT; USE NUM_INT
      IMPLICIT NONE
      DOUBLE PRECISION :: DIA(NNOD), P,Q,D,D1,D2,AR
      DOUBLE PRECISION :: sx,sy,syx,sz,pr
      DOUBLE PRECISION :: XG(*), YG(*), DISP(*), S1,S2, THETA
      DOUBLE PRECISION :: ST(4), E1(4), E2(4), E3(4), E4(4), E5(4)
      DOUBLE PRECISION :: P1(NNOD),P2(NNOD),P3(NNOD),P4(NNOD),P5(NNOD) 
      INTEGER :: I, IEL,ico(5,*), LJ(4), IX(*)
	   
C** FORM GLOBAL STIFFNESS MATRIX

      DIA(1:NNOD) = 0.D0
      P1(1:NNOD) = 0.D0;  P2(1:NNOD) = 0.D0
      P3(1:NNOD) = 0.D0;  P4(1:NNOD) = 0.D0;  P5(1:NNOD) = 0.D0
	
      DO IEL = 1,NEL
        LJ(1:4) = ICO(1:4,IEL) 
        AR = AREA(IEL)/4.D0
        ST(1:4) = AR;           CALL MAPLD(DIA,ST,LJ,4)
        E1(1:4) = AR*SXX(IEL);  CALL MAPLD(P1,E1,LJ,4)
        E2(1:4) = AR*SYY(IEL);  CALL MAPLD(P2,E2,LJ,4)
        E3(1:4) = AR*SXY(IEL);  CALL MAPLD(P3,E3,LJ,4)
        E4(1:4) = AR*SZZ(IEL);  CALL MAPLD(P4,E4,LJ,4)
        E5(1:4) = AR*PRESSURE(IEL); CALL MAPLD(P5,E5,LJ,4)
        
        P = (SXX(IEL)+SYY(IEL))/2.D0; 
        Q = SQRT((SXX(IEL)-SYY(IEL))**2/4.D0+SXY(IEL)**2)
        S1 = P + Q;  S2 = P - Q; D = SXX(IEL) - SYY(IEL) + 1.D-14
        THETA = (SXX(IEL)+SYY(IEL)+SZZ(IEL))/3.D0;
          
        WRITE(8,'(4I5,10E12.3)')ICO(1:4,IEL),SXX(IEL),SYY(IEL),SXY(IEL),
     +                           SZZ(IEL),PRESSURE(IEL),P,Q,S1,S2,THETA
      ENDDO

      REWIND 7
      DO I = 1, NNOD
	  IF(ABS(DIA(I)) < 1.D-08) DIA(I) = 1.0
        P1(I)=P1(I)/DIA(I);    P2(I)=P2(I)/DIA(I)
        P3(I)=P3(I)/DIA(I);    P4(I)=P4(I)/DIA(I)
        P5(I)=P5(I)/DIA(I)
    
        P = (P1(I)+P2(I))/2.D0; 
        Q = SQRT((P1(I)-P2(I))**2/4.D0+P3(I)**2)
        S1 = P + Q;  S2 = P - Q 
        D = P1(I) - P2(I) + 1.D-14
        THETA = 0
        D1 = 0.0; IF(IX(2*I-1).NE.0) D1 = DISP(IX(2*I-1))   
        D2 = 0.0; IF(IX(2*I).NE.0)   D2 = DISP(IX(2*I))
        
        WRITE(7,'(14E12.3)') XG(I),YG(I),D1,D2,P1(I),P2(I),P3(I),
     +                       P4(I),P5(I),P,Q,D1,D2,SQRT(D1**2+D2**2)
          
      ENDDO

      DO IEL = 1,NEL
        LJ(1:4) = ICO(1:4,IEL) 
        sx  = (p1(lj(1))+p1(lj(2))+p1(lj(3))+p1(lj(4)))/4.d0        
        sy  = (p2(lj(1))+p2(lj(2))+p2(lj(3))+p2(lj(4)))/4.d0        
        syx = (p3(lj(1))+p3(lj(2))+p3(lj(3))+p3(lj(4)))/4.d0        
        sz  = (p4(lj(1))+p4(lj(2))+p4(lj(3))+p4(lj(4)))/4.d0 
        pr  = (p5(lj(1))+p5(lj(2))+p5(lj(3))+p5(lj(4)))/4.d0       
        P   = (SX+SY)/2.D0; 
        Q   = SQRT((SX-SY)**2/4.D0+SYX**2)
        S1 = P + Q;  S2 = P - Q
        THETA = (SX+SY+SZ)/3.D0;
          
        WRITE(9,'(10E12.3)') SX,SY,SYX,SZ,pr,P,Q,s1,s2,theta
      ENDDO

      END SUBROUTINE SMOOTH	
		 
C...............................................................................

      SUBROUTINE STIFF_GRAVITY(IS,X,Y,S,ELOAD)
      USE MPROPERTY; USE NUM_INT; USE MPOINT      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),XI(3),W(3),H(7,7),Z(3,7),CMAT(3,3),C(3,7)
      DIMENSION H1(7,7),S(8,8),B(3,8),R(7,8),R1(7,8), ELOAD(10)
	  
      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/

      H(1:7,1:7)  = 0.D0;   H1(1:7,1:7)  = 0.D0
      R(1:7,1:8)  = 0.D0;   R1(1:7,1:8)  = 0.D0
      S(1:8,1:8)  = 0.D0;   ELOAD(1:8)   = 0.D0

C     THIS SUBROUTINE WILL CALCULATE THE MATRICES FOR EACH ELEMENT.  

      DO 26 I=1,3
      DO 27 J=1,3
        CALL SHAPE(XI(I),XI(J))
        CALL SHAPEST(X,Y,XL,YL,Z)    ! Z = P in Matlab code
        CALL BMATRX(X,Y,B,DET)
        
        CALL CMATRX(GMOD(IS),ANV(IS),CMAT)									   
        CALL MULT1(CMAT,Z,H1,C,3,7,3)
        DO 3 K=1,7
          DO 4 L=1,7
    4     H(K,L)=H(K,L)+(W(I)*W(J)*H1(K,L)*DET)
    3   CONTINUE

        DO K=1,7
          DO L = 1,8
            R1(K,L)  = DOT_PRODUCT(Z(:,K),B(:,L))  
          ENDDO
        ENDDO

        DO K=1,7
          DO L=1,8
            R(K,L)=R(K,L)+(W(I)*W(J)*R1(K,L)*DET)  ! R = Ls in Matlab code
          ENDDO
        ENDDO

   27 CONTINUE
   26 CONTINUE

      CALL INVERT (H,7)
      R1(1:7,1:8)  = 0.D0
      CALL MULT1(H,R,S,R1,7,8,7)     ! R(transpose)*H*R

      RETURN
      END SUBROUTINE STIFF_GRAVITY

C...............................................................................

      SUBROUTINE STIFFJ(IEL,X,Y,S)
      USE mproperty; USE mpoint; USE joint
      IMPLICIT None
      DOUBLE PRECISION :: X(*),Y(*),S(4,4),D(2,2)
      DOUBLE PRECISION :: B(2,4),C(2,4), DET
      INTEGER :: IEL

C     THIS SUBROUTINE WILL CALCULATE THE STIFFNESS MATRIX FOR EACH JOINT ELEMENT.  

      S(1:4,1:4) = 0.d0
      D(1,1) =  KNS(IEL);  D(1,2) = 0.d0
      D(2,1) =  0.d0;      D(2,2) = KNN(IEL)
      CALL BMATRXJ(X,Y,B,DET)
      CALL MULT1(D,B,S,C,2,4,2)
      S(1:4,1:4) = S(1:4,1:4)*DET

      RETURN
      END SUBROUTINE STIFFJ

C..........................................................................

      SUBROUTINE UPDATE(XG,YG,ICO,IX,DISP)
      USE GCONTROL; USE MPROPERTY; USE NUM_INT; USE MPOINT; USE JOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: X(8),Y(8),XG(*),YG(*),V(10),V1(8)
      DOUBLE PRECISION :: DISP(*),ST(4)
      DOUBLE PRECISION :: B(7),C(7,8)
      DOUBLE PRECISION :: XL,YL,DET,US(4),ENN,GNS
      DOUBLE PRECISION :: H(7,7),Z(3,7),CMAT(3,3),Z1(3,7),R1(7,8)
      DOUBLE PRECISION :: BMAT(3,8),R(7,8),XI(3),W(3),H1(7,7),BJ(2,4)
      INTEGER :: I, IEL,J,K, L, IS,ico(5,1)
      INTEGER :: IX(*),LJ(18)
      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/

C**   UPDATE STRESSES

      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC(IEL),YMC(IEL))

        V(1:10) = 0.D0
        DO J=1,10
          IF(LJ(J) > 0) V(J)=DISP(LJ(J))
        ENDDO
        V1(1:8) = V(1:8)

        H(1:7,1:7)   = 0.D0;  H1(1:7,1:7)  = 0.D0
        R(1:7,1:8)   = 0.D0;  R1(1:7,1:8)  = 0.D0


        DO I=1,3
        DO J=1,3
          CALL SHAPE(XI(I),XI(J))
          CALL SHAPEST(X,Y,XL,YL,Z)
          CALL BMATRX(X,Y,BMAT,DET)
          CALL CMATRX(GMOD(IS),ANV(IS),CMAT)
          CALL MULT1(CMAT,Z,H1,Z1,3,7,3)

          DO 3 K=1,7
            DO 4 L=1,7
    4        H(K,L)=H(K,L)+(W(I)*W(J)*H1(K,L)*DET)
    3     CONTINUE

          DO K=1,7
            DO L = 1,8
              R1(K,L)  = DOT_PRODUCT(Z(:,K),BMAT(:,L))
            ENDDO
          ENDDO

          DO K=1,7
            DO L=1,8
              R(K,L)=R(K,L)+(W(I)*W(J)*R1(K,L)*DET)
            ENDDO
          ENDDO

        ENDDO
      ENDDO

      CALL INVERT (H,7)

        DO K=1,7
          DO L = 1,8
            C(K,L)  = DOT_PRODUCT(H(K,:),R(:,L))
          ENDDO
        ENDDO

        DO K=1,7
            B(K)  = DOT_PRODUCT(C(K,:),V1(:))
        ENDDO

        ST(1) = B(1);      
        ST(2) = B(4);     
        ST(3) = B(7);         
        ST(4) = ANV(IS)*(ST(1)+ST(2));   

        CALL GLOLOC(SXX,SYY,SZZ,SXY,ST,IEL,0)

      ENDDO

      IF(NELJ > 0) then
        DO IEL=1,NELJ
          CALL LOCALJ(XG,YG,X,Y,IEL,LJ,IX)         
          us(1:4) = 0.d0
          DO J=1,4
            IF(LJ(J) > 0)  US(J) = disp(LJ(J))
          Enddo

          call bmatrxJ(X,Y,BJ,DET)

C** CALCULATE THE INTERFACE STRAIN 

          GNS  = Dot_Product(BJ(1,:),US(:))
          ENN  = Dot_Product(BJ(2,:),US(:))

          TNS(IEL) = TNS(IEL) + KNS(IEL)*GNS
          TNN(IEL) = TNN(IEL) + KNN(IEL)*ENN
        ENDDO  
      ENDIF
	  
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
		  
C......................................................................

      SUBROUTINE BANDWH(ICO,JX,LJ)
      USE gcontrol
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ICO(5,*),JX(*),LJ(*)
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
      END


C.........................................................................

      SUBROUTINE INVERT (A,N)
C
C*****************************************
C     COMPUTES THE INVERSE MATRIX OF A
C*****************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),K(25)
C
      DET=1.D0
      DO 5 I=1,N
   5  K(I)=I
      DO 80 II=1,N
      DO 10 I=II,N
      PIV=A(I,II)
      IF (DABS(PIV).GT.1.D-12) GOTO 20
  10  CONTINUE
      WRITE(6,*)'MATRIX IS SINGULAR'
      STOP
  20  DET=DET*PIV
      IF(I.EQ.II) GOTO 40
      I1=K(II)
      K(II)=K(I)
      K(I)=I1
      DO 30 J=1,N
      C=A(I,J)
      A(I,J)=A(II,J)
  30  A(II,J)=C
      DET=-DET
  40  C=1.D0/PIV
      A(II,II)=1.D0
      DO 50 J=1,N
  50  A(II,J)=A(II,J)*C
      DO 70 I=1,N
      IF(I.EQ.II) GOTO 70
      C=A(I,II)
      A(I,II)=.0D0
      DO 60 J=1,N
  60  A(I,J)=A(I,J)-C*A(II,J)
  70  CONTINUE
  80  CONTINUE
      DO 120 J=1,N
      DO 90 J1=J,N
      JJ=K(J1)
      IF(JJ.EQ.J) GOTO 100
  90  CONTINUE
 100  IF (J.EQ.J1) GOTO 120
      K(J1)=K(J)
      DO 110 I=1,N
      C=A(I,J)
      A(I,J)=A(I,J1)
 110  A(I,J1)=C
 120  CONTINUE
      RETURN
      END SUBROUTINE INVERT

C...........................................................................

      SUBROUTINE LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX,XMC,YMC)
      IMPLICIT NONE
      DOUBLE PRECISION :: XG(*),YG(*),X(*),Y(*),XMC,YMC
      INTEGER ::  IS, IEL, I, J, J1, J2, NNODEL, NVAR, ICOO 
      INTEGER ::  LJ(*),IX(*),ico(5,*)

      NVAR = 2
      NNODEL = 4

      DO I=1,NNODEL
        ICOO=ICO(I,IEL)
        X(I)=XG(ICOO)-XMC
        Y(I)=YG(ICOO)-YMC
      ENDDO

      IS=ICO(NNODEL+1,IEL)
      DO J=1,NNODEL
        J1=(J-1)*NVAR
        J2=NVAR*(ICO(J,IEL)-1)
        DO I=1,NVAR
           LJ(I+J1)=IX(J2+I)
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE LOCAL

C...........................................................................

      SUBROUTINE LOCALB(XG,YG,X,Y,IEL,LJ,IX)
      USE TRACTIONS
      IMPLICIT NONE
      DOUBLE PRECISION XG(*),YG(*),X(*),Y(*)
      INTEGER :: IEL, I, J, J1, J2, NVAR, NNODEL, LJ(*),IX(*), IC(2)
   
      NVAR = 2
      NNODEL = 2

      DO I=1,NNODEL
       IC(I)=ICOB(I,IEL)
       X(I)=XG(ICOB(I,IEL))
       Y(I)=YG(ICOB(I,IEL))
      ENDDO

      DO J=1,NNODEL
         J1=NVAR*(J-1)
         J2=NVAR*(IC(J)-1)
        DO I=1,NVAR
          LJ(I+J1)=IX(J2+I)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE LOCALB

C...........................................................................

      SUBROUTINE LOCALJ(XG,YG,X,Y,IEL,LJ,IX)
      USE joint
      Implicit None
      DOUBLE PRECISION XG(*),YG(*),X(*),Y(*)
      Integer :: iel, i, j, j1, j2, nvar, nnodel, LJ(*),IX(*)
   
      NVAR = 2
      nnodel = 2
      X(1:NNODEL)=XG(ICOJ(1:NNODEL,IEL))
      Y(1:NNODEL)=YG(ICOJ(1:NNODEL,IEL))

      DO J=1,nnodel
         J1=NVAR*(J-1)
         J2=NVAR*(ICOJ(J,IEL)-1)
        DO I=1,NVAR
          LJ(I+J1)=IX(J2+I)
        enddo
      enddo
      RETURN
      END SUBROUTINE LOCALJ

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

      SUBROUTINE MAPST(A,S,NVEL,LJ,LBAND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),S(NVEL,NVEL),LJ(*)
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
   11 CONTINUE
   12 CONTINUE
      RETURN
      END

C...............................................................................

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
      END

C..................................................................

      SUBROUTINE GAUSS
      USE NUM_INT
      DOUBLE PRECISION :: A 

      A = 1.D0/SQRT(3.D0)
C      A = 0.5D0
  
      XS(1) = -A;    XS(2) =  A ;   XS(3) = A;    XS(4) = -A
      XT(1) = -A;    XT(2) = -A;    XT(3)  = A;   XT(4) = A
      WF(1) = 1.D0;  WF(2) = 1.D0;  WF(3) = 1.D0; WF(4) = 1.D0

      RETURN
      END


C..................................................................
