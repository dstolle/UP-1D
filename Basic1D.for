      MODULE GCONTROL
       INTEGER, PARAMETER :: NNG = 1000, NELT  = 1000
       INTEGER, PARAMETER :: NFD = 1000, MSTIF = 100000 
       INTEGER :: NEL,NNOD,NVAR,NNODEL,NVEL,NMAT,NNET,LBAND,NVA
       INTEGER :: IPS,MATR,NCOD,NCON,NLN,IGR,NSTEP,NITER
       INTEGER :: NVAp,NNETp,LBANDp,idisp = 10, NPT =100
       DOUBLE PRECISION :: Tr, Cvmax, time, dtime ! time step, Tr = applied pressure
       DOUBLE PRECISION :: DH,Rmin, Hmin	! Rmin = location of zero excess pressure
                                          ! Hmin = smallest element length
                                          ! DH = Height of element 
       CHARACTER (LEN = 5) :: INFILE
      END MODULE GCONTROL

      MODULE MPROPERTY
        DOUBLE PRECISION, DIMENSION ::GMOD,EMOD,ANV,Em
        DOUBLE PRECISION, DIMENSION :: GRR, Kw,Cv,BULK
        DOUBLE PRECISION :: SIG0, GRw = 10, beta
      END MODULE MPROPERTY

      MODULE MPOINT
        INTEGER, PARAMETER :: IPNT=1000
        DOUBLE PRECISION, DIMENSION (IPNT)::GAMA,UW
        DOUBLE PRECISION, DIMENSION (IPNT)::SXX,SYY,SXY,SZZ
      END MODULE MPOINT


C******************************************************************

      PROGRAM Basic1D
      USE GCONTROL; USE MPROPERTY; USE NUM_INT;USE TRACTIONS; USE MPOINT
				
C       ************************************************************
C       *                                                          *
C       *      2-NODED STRESS ANALYSIS FINITE ELEMENT PROGRAM      *
C       *        CONSOLIDATION WITH P CONSTANT and U linear        *
C       *                    This is linear elastic                *
C       *                       (February 2022)                    *
C       ************************************************************

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION (NFD) :: TLOAD,PLOAD,TDISP !,BB
      DOUBLE PRECISION :: XG(NNG),YG(NNG),X(4),Y(4)
      DOUBLE PRECISION :: GSTIF(MSTIF),CC(NNG),CN(8)
      DOUBLE PRECISION :: DET, RATIO
      DOUBLE PRECISION :: P0(NNG),Ps(NELT),PT(NELT),dPT(NELT)
      INTEGER :: I,IEL,NCN,istep,IS
      INTEGER :: ICO(3,NELT),IX(NFD),IXp(NFD),LJ(4)

      WRITE(*,'(A,\)') ' INPUT FILE NAME = '
      READ(*,'(A)') INFILE

      OPEN(2,FILE=INFILE//'.DIS',FORM='FORMATTED')  ! Initial Displacements
      OPEN(4,FILE=INFILE//'.STR',FORM='FORMATTED')  ! Initial stresses
      OPEN(5,FILE=INFILE//'.DAT',FORM='FORMATTED')  ! Basic FEM input
      OPEN(6,FILE=INFILE//'.OUT',FORM='FORMATTED')  ! Output at end of analysis
      OPEN(7,FILE=INFILE//'.NND',FORM='FORMATTED')  ! For Cplot Input
      OPEN(8,FILE=INFILE//'.IEL',FORM='FORMATTED')  ! For Cplot Input
 
C** INPUT THE INFORMATION NEEDED FOR THE REQUESTED FINITE ELEMENT ANALYSIS

      CALL INPUT(XG,ICO,IX,IXp)

C*******************************************************************

C**  STEP 1: INITIALIZES ARRAYS INCLUDING INITIAL STRESSES

      dPT(1:NEL) = 0.D0;  P0(1:NNOD) = 0.D0;

C**  INCLUDES SETUP OF MATRICES FOR MECHANICAL ANALYSES
C**   Form load vector

      TLOAD(1) = 	TR*R*DH               ! Surface load
      TDISP(1:NNET) = 0.D0;  TLOAD(1:NNET) = 0.D0
 

      CALL MATRX(GSTIF,CC,ICO,XG,DH,IX) ! Set up matrices





C**   Initial effective stresses at centroid, plus pressures
      DO I =1, NEL
        READ(4,*) SXX(I),SYY(I),SXY(I),SZZ(I),UW(I)	     
        Ps(I)  = -(SXX(I)+SYY(I)+SZZ(I))/3.d0	! Effective pressure p'
        PT(I)  = Ps(I)+ UW(I)  ! Compression positive, total pressure P     
      ENDDO                    ! Tension +ve for Sij    

      CALL Pintegral(P0,UW,ICO,XG,YG,IX,1) ! nodal pressure given UW

      DO I = 1, NNOD		
        if(IXp(I)== 0) P0(I) = 0.D0         ! Essential BC for P0
      ENDDO
   
      WRITE(6,101)
      CALL OUTPUT(XG,YG,ICO,TDISP,IX,P0)

C*******************************************************************

C**  STEP 2: SETUP MATRICES FOR MECHANICAL & SEEPAGE ANALYSES
C**   Form load vector

      CALL LOAD(TLOAD,XG,YG,IX)         ! Surface load
      CALL MATRX(GSTIF,CC,ICO,XG,YG,IX) ! Set up matrices

C**    CHECK F - Ka calculation
c      CALL MatrixXa(BB,TDISP,GSTIF,NNET,LBAND) 
c      BB(1:nnet) = Tload(1:nnet) - BB(1:nnet)      

C**   Determine the time step

      dtime = Hmin**2/(4*Cvmax)		 ! Critical time step

      WRITE(*,'(A\,E12.5,/)')   ' Critical time step  =  ', dtime
      WRITE(*,'(A\)') ' Provide dtime     =  ' 
      READ(*,*) dtime           
      WRITE(*,'(A\)') ' Provide nstep     =  '
      READ(*,*) nstep

      WRITE(6,'(/,A\,E12.5)')   ' Time step  =  ', dtime      

      WRITE(6,'(A\,I7,/)')    ' Number of steps   =  ', nstep
      WRITE(6,'(A\,E12.5,//)')' Maximum time      =  ', nstep*dtime

      WRITE(*,'(//)')                   

C*******************************************************************

C**  STEP 3:  LOOP FOR TIME STEPPING

      time = 0.d0
      DO istep = 1, nstep
        time = time + dtime

C**   UPDATE PORE PRESSURES
        
        CALL SEEPAGE(CC,P0,dPT,ICO,IX,IXp,XG,YG) 
	                                             	  	               
C**   MECHANICAL STEP BASED ON EFFECTIVE STRESS CHANGES

        PLOAD(1:NNET) = 0.D0
        DO IEL = 1, NEL       ! Qt*uw 
          CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX)
          CALL PVEC(X,Y,CN)
          CN(1:8) =  UW(IEL)*CN(1:8)
          CALL MAPLD(PLOAD,CN,LJ,NVEL)
        ENDDO
    
C**   Displacements associated with effective stresses
C**   Calculate net load associated with displacement calculation

        TDISP(1:nnet) = Tload(1:nnet) + PLOAD(1:nnet)

        DET=1.D-8  
        RATIO=1.D-8  
        CALL CHOLSQ(GSTIF,TDISP,NNET,LBAND,2,RATIO,DET,NCN)
        CALL UPDATE(XG,YG,ICO,IX,TDISP,Ps,dPT)
        IF((istep/NPT)*NPT.EQ.istep) write(*,'(I10,3E12.5)') istep,time,
     +                                                UW(1),TDISP(idisp)
      ENDDO

      CALL OUTPUT(XG,YG,ICO,TDISP,IX,P0)
      CALL SMOOTH(ICO,XG,YG,TDISP,IX)

      WRITE(*,*) '*** PRESS ANY KEY ***'
      READ(*,*) 					
   
      WRITE(6,102)

      CLOSE(2);  CLOSE(4);  
      CLOSE(5);  CLOSE(6);  CLOSE(7);  CLOSE(8);

      STOP

  101  FORMAT(//,' INITIAL MATERIAL POINT DISPLACEMENTS & STRESSES',
     +        /,' ************************************************')
  102 FORMAT(///12X,'*************** END OF ANALYSIS ***************'/)
 104  FORMAT(I5,6X,6E13.5)
      END

C.............................................................

      SUBROUTINE BMATRX(X,Y,B,DET,R)
      USE gcontrol; USE num_int
      Implicit None
      DOUBLE PRECISION :: X(4),Y(4),B(4,*),R
      DOUBLE PRECISION :: J11,J12,J21,J22, DET, AI(2,2)

      B(1:4,1:8) = 0.d0
      J11 = Dot_Product(ans,x); J12 = Dot_Product(ans,y)
      J21 = Dot_Product(ant,x); J22 = Dot_Product(ant,y)

      DET = J11*J22 - J12*J21 
      AI(1,1) =  J22/DET;  AI(1,2) = -J12/DET                   
      AI(2,1) = -J21/DET;  AI(2,2) =  J11/DET                   
      
      B(1,1:7:2)=AI(1,1)*ANS(1:4)+AI(1,2)*ANT(1:4)       
      B(2,2:8:2)=AI(2,1)*ANS(1:4)+AI(2,2)*ANT(1:4)       
      B(3,1:7:2)=B(2,2:8:2)  
      B(3,2:8:2)=B(1,1:7:2)

      if(igrad > 0) then		! SELECTIVE INTEGRATION FOR SHEAR
        CALL AVGRAD(X,Y) 
        B(3,1:7:2) = YN(1:4)
        B(3,2:8:2) = XN(1:4)
      ENDIF

      R = 1.d0
      if (IPS > 0) then
        R =  Dot_Product(an,x)
        B(4,1:7:2)=AN(1:4)/R
      endif

      RETURN
      END SUBROUTINE BMATRX 

C.........................................................................
		   
      SUBROUTINE CONDUCTANCE(IS,X,Y,H)
      USE GCONTROL; USE MPROPERTY; USE MPOINT;; USE num_int 
      IMPLICIT NONE
      DOUBLE PRECISION :: H(4,4),X(4),Y(4),W(3),XI(3)
      DOUBLE PRECISION :: R,DET,B(2,4),J11,J12,J21,J22,AI(2,2)
      INTEGER :: I, J, IS, K ,L

      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/

C** FORM LOCAL CONDUCTANCE MATRIX

      H(1:4,1:4) =0.D0

      DO I=1,3		  ! INTEGRATE OVER INTEGRATION POINTS
        DO J=1,3
           CALL SHAPE(XI(I),XI(J))
           B(1:2,1:4) = 0.D0      ! FORM B-MATRIX FOR SEEPAGE
           J11 = DOT_PRODUCT(ANS,X); J12 = DOT_PRODUCT(ANS,Y)
           J21 = DOT_PRODUCT(ANT,X); J22 = DOT_PRODUCT(ANT,Y)
           DET = J11*J22 - J12*J21 
           AI(1,1) =  J22/DET;  AI(1,2) = -J12/DET                   
           AI(2,1) = -J21/DET;  AI(2,2) =  J11/DET                    
           B(1,1:4)=AI(1,1)*ANS(1:4)+AI(1,2)*ANT(1:4)       
           B(2,1:4)=AI(2,1)*ANS(1:4)+AI(2,2)*ANT(1:4)

           R = 1.d0
           if (IPS > 0) R =  Dot_Product(an,x)

           DET = Cv(IS)*W(I)*W(J)*DET*R 	  
	 
	     DO K = 1, 4		!	Bt*B
             DO L = 1, 4 
               H(K,L) = H(K,L) + (B(1,K)*B(1,L)+B(2,K)*B(2,L))*DET
             ENDDO
           ENDDO
        ENDDO    !J
      ENDDO    !I

      RETURN
      END SUBROUTINE CONDUCTANCE		   

C............................................................

      SUBROUTINE DISPL(DIS,IX,XG,YG,P1)
      USE GCONTROL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DIS(*),IX(*),D(4),XG(*),YG(*),P1(*)

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

      WRITE(6,103) I,XG(I),YG(I),(D(K),K=1,NVAR),P1(I)

 1    CONTINUE
      RETURN
  102 FORMAT(//,1X,'N O D A L    D I S P L A C E M E N T S.........',
     +  //,1X,'NODE   X-COORDINATE   Y-COORDINATE',
     +     1X,'DISPLACEMENT-U  PRESSURE',/)
  103 FORMAT(I5,5E15.6)
      END SUBROUTINE DISPL 

C........................................................................

      SUBROUTINE DMATRX(E,ANU,D)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION D(4,4)
      D(1:4,1:4) = 0.D0

        E1=E/(1.D0-ANU*ANU)
        ANV=ANU/(1.D0-ANU)
        D(1,1)=E1/(1.D0-ANV*ANV); D(1,2)=ANV*D(1,1); D(1,4)=ANV*D(1,1)  
        D(2,1)=D(1,2);            D(2,2)=D(1,1);	 D(2,4)=D(1,2)
        D(3,3)=E/(2.D0*(1.D0+ANU))
        D(4,1)=D(1,4);            D(4,2)=D(2,4);	 D(4,4)=D(1,1)

      RETURN
      END SUBROUTINE DMATRX

C............................................................

      SUBROUTINE INPUT(XG,YG,ICO,IX,IXp)
      USE GCONTROL;USE MPROPERTY;USE NUM_INT;USE TRACTIONS;  
      USE MPOINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XG(*),YG(*),ico(2,*),IX(*),IXp(*),LJ(10)
      CHARACTER *64 TITLE

C** READ IN THE MASTER CONTROL INFORMATION

      READ(5,'(A)')     TITLE
      READ(5,*)  NEL,NNOD,Rmin,Rmax,DH,Tr
C
**   IPS 0 > 1: axisymmetric

      WRITE(6,100)

      CALL LAYOUT(XG,YG,ICO,IX,IXp)
      CALL BANDWH(ICO,IX,LJ)
      NVA = (LBAND+1)*NNET 

      WRITE(6,'(1X,A)') TITLE

      WRITE(6,102) NEL,NNOD,Rmin,Rmax,DH,Tr    

      WRITE(6,110) NVEL,NNET,LBAND,NVA

C**  MATERIAL PROPERY SETS

      Cvmax = 1.D-10
      WRITE(6,115)
      READ(5,*) EMOD,ANV,GRR,Kw
      GMOD = EMOD/(1.D0+ANV)/2.D0
      Em   = 2.D0*GMOD*(1- ANV)/(1.D0-2.D0*ANV)
      Bulk = EMOD/(3.d0*(1.D0-2.D0*ANV))
      Cv	  = Em*Kw/GRw
      if(Cv > Cvmax) Cvmax = Cv
      WRITE(6,116) EMOD,ANV,GRR,Em,Bulk,Kw,Cv
      ENDDO   

  100 FORMAT(/,5x,'***** UP-CONSOLIDATION: AXI-SYMMETRIC *****',/)    
  102 FORMAT( /,5X,'TOTAL NO. OF ELEMENTS             NEL   = ',I5,
     +        /,5X,'TOTAL NUMBER OF ELEMENTS         NNOD   = ',I5,
     +        /,5X,'BOREHOLE RADIUS (RMIN)..................= ',E12.5,       
     +        /,5X,'DOMAIN RADIUS (RMAX)....................= ',E12.5,       
     +        /,5X,'DOMAIN THICKNESS (DH)...................= ',E12.5,       
     +        /,5X,'BOREHOLE WALL PRESSURE (Tr).............= ',E12.5)

  110 FORMAT(/,5X,'NO. OF VARIABLES PER ELEMENT      NVEL  =',I5,
     +       /,5X,'NET DEGREES OF FREEDOM            NDEG  =',I5,
     +       /,5X,'HALF BANDWIDTH EXCLUDING DIAG.    LBAND =',I5,
     +       /,5X,'SIZE OF GLOBAL STIFFNESS MATRIX   NVA   =',I8,/)
 115  FORMAT(/,1X,'M A T E R I A L    P R O P E R T I E S ')
 116  FORMAT(//,'PROPERTIES FOR MATERIAL SET ',
     +        /,'*************************** ',
     +        /,'ELASTIC MODULUS (EMOD)............= ',E12.5,
     +        /,'POISSON RATIO (ANV)...............= ',E12.5,
     +        /,'UNIT WEIGHT (GRR)................ = ',E12.5,
     +        /,'CONSTRAINED MODULUS (Em)......... = ',E12.5,
     +        /,'BULK MODULUS (BULK).............. = ',E12.5,
     +        /,'PERMEABILITY (Kw)................ = ',E12.5,
     +        /,'COEFFICIENT OF CONSOLIDATION (Cv) = ',E12.5,/)
    2 FORMAT(I5,5X,2I5,5X,6E10.3)
      END SUBROUTINE INPUT

C.........................................................................

      SUBROUTINE LAYOUT(XG,ICO,IX,IXp)
      USE GCONTROL; 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XG(*),IX(*),IXp(*),ico(2,*)

      NNODEL = 2;  NNN=NNODEL+1; NVAR = 1;  NVEL=NVAR*NNODEL;	  

      WRITE(6,41) NEL,NNOD,NVAR,NNODEL

      DL = (Rmax - Rmin)/NEL
      HMIN = DL

      DO I = 1, NNOD	         ! create nodes
        XG(I) = (I-1)*DL + Rmin;  IX(I) = 1              
      ENDDO    

      DO IEL = 1, NEL
        ICO(1,IEL) = IEL; ICO(2,IEL) = IEL + 1        
      ENDDO        

      NMAT=NVAR*NNOD
    
      IXp(1:NNOD) = 1
      IXp(1) = 0  ! ZERO EXCESS PORE PRESSURE AT FACE   

C     WRITE KEY INPUT

      WRITE(6,100)
      WRITE(6,42)
      DO I=1,NNOD
        WRITE(6,44) I,XG(I),IX(I),IXp(I)
      ENDDO
      
      WRITE(6,47)

      DO IEL = 1, NEL
        WRITE(6,48) IEL,(ICO(J,IEL),J=1,NNODEL)
      ENDDO

      NNET=0
      DO 12 I=1,NMAT
      IF(IX(I).GT.0) THEN
        NNET=NNET+1
        IX(I)=NNET
      ELSE
        IX(I)=0
      ENDIF
   12 CONTINUE

   41 FORMAT(//,' TOTAL NO. OF ELEMENTS..........',I5,/,
     +          ' NUMBER OF NODES................',I5,/,
     +          ' VARIABLES PER NODE.............',I5,/,
     +          ' NODES PER ELEMENT..............',I5,/)
   42 FORMAT(/,3X,' NODE',6X,' R-CORD',7X,' U   P')	
   44 FORMAT(1X,I5,5X,F10.4,5X,2I4)   
   47 FORMAT(//,5X,'ELEMENT',9X,'NODE NUMBERS',7X,'DRadius',//)
   48 FORMAT(5X,I5,6X,2I5,E12.4)
 100  FORMAT(//,' NODAL COORDINATES AND BOUNDARY CONDITIONS',
     +        /,' *****************************************')
      END SUBROUTINE LAYOUT
	  
C............................................................
			   
      SUBROUTINE MATRX(GSTIF,CC,ICO,XG,YG,IX)
      USE GCONTROL; USE MPROPERTY; USE MPOINT; USE NUM_INT;
      IMPLICIT NONE
      DOUBLE PRECISION :: GSTIF(*),ESTIF(8,8),CC(*),XG(*),YG(*)
      DOUBLE PRECISION :: ECC(4),F(nnet)
      DOUBLE PRECISION :: DET, RATIO, X(4),Y(4)
      INTEGER :: IEL, IS, ICO(5,*),IX(*),LJ(8),LJC(4),NCN

C** FORM GLOBAL STIFFNESS MATRIX
	
      GSTIF(1:NVA) = 0.D0;  CC(1:NNOD) = 0.D0;  F(1:NNET) = 0.D0

C**   Element  stiffness matrix H; Lumped Capacitance Matrix ECC 
      DO IEL=1,NEL
         CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX)
         CALL STIFF(IS,X,Y,ESTIF,ECC)			  ! 
         CALL MAPST(GSTIF,ESTIF,NVEL,LJ,LBAND)
         LJC(1:4) = ICO(1:4,IEL)
         CALL MAPLD(CC,ECC,LJC,4)	
      ENDDO

      DET=1.D-8
      RATIO=1.D-8
      CALL CHOLSQ(GSTIF,F,NNET,LBAND,1,RATIO,DET,NCN)
      IF(DET.LE.0.D0) WRITE(6,'(//,A\,E12.5)') ' DET < 0  =====>>> ',DET
	
      END SUBROUTINE MATRX		   
	
C............................................................

      SUBROUTINE OUTPUT(XG,YG,ICO,TDISP,IX,P1)
      USE GCONTROL; USE MPROPERTY; USE NUM_INT; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: X(4),Y(4),XG(*),YG(*),P1(*)
      DOUBLE PRECISION :: TDISP(*),DS(4),ST(4),XXX,YYY
      INTEGER :: IEL,IS
      INTEGER :: IX(*),ICO(5,*), LJ(12)

C** PRINT DISPLACEMENT

      CALL DISPL(TDISP,IX,XG,YG,P1)

      WRITE(6,101)
      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX)
        CALL GLOLOC(SXX,SYY,SZZ,SXY,DS,IEL,1)
        ST(1) = DS(1)- UW(IEL); ST(2) = DS(2)- UW(IEL);
        ST(3) = DS(3);          ST(4) = DS(4)- UW(IEL);
        XXX = (X(1)+X(2)+X(3)+X(4))/4.D0	
        YYY = (Y(1)+Y(2)+Y(3)+Y(4))/4.D0	 
        WRITE(6,104) IEL,XXX,YYY,ST(1:4),UW(IEL)
      ENDDO

  101 FORMAT(//,' MATERIAL POINT STRESSES',
     +        /,' ***********************',
     +   //,1X,'MPN    REFERENCE   LOCATION   STRESS -11   STRESS -22',
     +     '   STRESS -12   STRESS -33     PRESSURE',//)
  104 FORMAT(I5,3X,E9.3,2X,E9.3,7E13.4)
      END SUBROUTINE OUTPUT

C.........................................................................
                        	   
      SUBROUTINE Pintegral(P,Pc,ICO,XG,YG,IX,Nopt)
      USE GCONTROL; USE NUM_INT; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: PcL(4),XG(*),YG(*),P(*),Pc(*),DIA(NNOD)
      DOUBLE PRECISION :: R,DET,ECC(4),X(4),Y(4),W(3),XI(3)
      DOUBLE PRECISION :: J11, J12, J21, J22
      INTEGER :: I,J,IEL,IS,ICO(5,*),IX(*),LJ(8),LJC(4),Nopt

      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/
			 	
      DIA(1:NNOD) = 0.D0; P(1:NNOD)  = 0.D0  ! Initialize nodal P
	 
      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX)
        LJC(1:4) = ICO(1:4,IEL)
        ECC(1:4) = 0.D0;

        DO I=1,3
          DO J=1,3
            CALL SHAPE(XI(I),XI(J))
            J11 = Dot_Product(ans,x); J12 = Dot_Product(ans,y)
            J21 = Dot_Product(ant,x); J22 = Dot_Product(ant,y)
            DET = J11*J22 - J12*J21 

            R = 1.d0
            if (IPS > 0) R =  Dot_Product(an,x)
            DET = W(I)*W(J)*DET*R      
            ECC(1:4)   = ECC(1:4)  +  AN(1:4)*DET
          ENDDO    !J
        ENDDO    !I
                      
        PcL(1:4) = ECC(1:4)*Pc(IEL)					

        CALL MAPLD(DIA,ECC,LJC,4)			
        CALL MAPLD(P,PcL,LJC,4)			
      ENDDO
    
      IF(Nopt > 0) THEN
        DO I = 1, NNOD
          P(I) = P(I)/DIA(I)
        ENDDO 
      ENDIF

      RETURN        	
      END SUBROUTINE Pintegral

C............................................................
	
      SUBROUTINE PVEC(X,Y,CN)
      USE NUM_INT; USE MPOINT; USE num_int   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),CN(*),BMAT(4,8),W(3),XI(3)
      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/

C     THIS SUBROUTINE CALCULATES THE UW COEFFICIENTS: CN

      CN(1:8) = 0.D0;
	 
      DO I=1,3
        DO J=1,3
           CALL SHAPE(XI(I),XI(J))
           CALL BMATRX(X,Y,BMAT,DET,R)
           DET = W(I)*W(J)*DET*R      
           CN(1:7:2)  = CN(1:7:2) + (BMAT(1,1:7:2)+BMAT(4,1:7:2))*DET
           CN(2:8:2)  = CN(2:8:2) + (BMAT(2,2:8:2))*DET
        ENDDO    !J
      ENDDO    !I

      RETURN
      END SUBROUTINE PVEC     		   

C............................................................
	      
      SUBROUTINE SEEPAGE(CC,P0,dPT,ICO,IX,IXp,XG,YG)
      USE GCONTROL; USE MPROPERTY; USE MPOINT; USE NUM_INT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ::  RLOAD(NNOD),Fi(NNOD),dPT(*),R(NNOD) 
      DOUBLE PRECISION ::  CC(*),X(4),Y(4),P0(*),XG(*),YG(*)                    
      DOUBLE PRECISION ::  FL(4),P(4),dP(4),H(4,4),P1(nnod)
      INTEGER:: ICO(5,*), I, K, IEL, LJ(4), LJC(4), IX(*), IXp(*) 

      RLOAD(1:NNOD) = 0.d0; Fi(1:NNOD) = 0.d0

C**   UPDATE PORE PRESSURE

C**   Force due to change in total pressure: Step 1

      CALL Pintegral(RLOAD,dPT,ICO,XG,YG,IX,0) ! RLo = A*dP
      R(1:NNOD)= RLOAD(1:NNOD)

C**   Internal force due to pore pressure: Step 2	 

      DO IEL = 1, NEL	        
        LJC(1:4) = ICO(1:4,IEL); IS = ICO(5,IEL)
        X(1:4)=XG(LJC(1:4));  Y(1:4)=YG(LJC(1:4))        
        CALL CONDUCTANCE(IS,X,Y,H)
        P(1:4)= P0(LJC(1:4))   
        DO K=1,4
          FL(K) =  DOT_PRODUCT(H(K,:),P(:))	  !  H*p0
        ENDDO
        CALL MAPLD(Fi,FL,LJC,4)         
      ENDDO

C**   Unbalanced load, RL = RLo - dt*H*p0

      RLOAD(1:NNOD)=RLOAD(1:NNOD)-dtime*Fi(1:NNOD)

      DO I = 1, NNOD
        IF(IXp(I) > 0) then
          P1(I) = P0(I) + RLOAD(I)/CC(I)
        ELSE
          P1(I) = 0.d0
        ENDIF

      ENDDO

C**   UPDATE ELEMENT PORE PRESSURE

      DO IEL = 1, NEL
        LJ(1:4) = ICO(1:4,IEL)
        X(1:4)=XG(LJ(1:4)); Y(1:4)=YG(LJ(1:4))                
        P(1:4)   =  P1(LJ(1:4))
        dP(1:4)  =  P1(LJ(1:4)) - P0(LJ(1:4))
        UW(IEL)  = (P(1)+P(2)+P(3)+P(4))/4.D0 ! UW at centroid      
C**   Stage 1 update for dPT
        dPT(IEL) = (dP(1)+dP(2)+dP(3)+dP(4))/4.D0 
      ENDDO                                           

      P0(1:NNOD) = P1(1:NNOD)                !  update p0
      
      END SUBROUTINE SEEPAGE

C............................................................

      SUBROUTINE STIFF(IS,X,Y,S,ECC)
      USE MPROPERTY; USE NUM_INT; USE MPOINT; USE num_int   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),S(8,8),BDB(8,8),DMAT(4,4),ECC(*)
      DIMENSION BMAT(4,8),C(4,8),W(3),XI(3),CN(8)
      DATA  W/ 0.55555555555556,0.88888888888889,0.55555555555556/
      DATA XI/-0.77459666924148,0.00000000000000,0.77459666924148/

C     THIS SUBROUTINE WILL CALCULATE THE STIFFNESS MATRIX FOR EACH ELEMENT.  

      AR = 0.D0;  CN(1:8) = 0.D0; S(1:8,1:8) = 0.D0; ECC(1:4) = 0.D0

      CALL DMATRX(EMOD(IS),ANV(IS),DMAT)
	 
      DO I=1,3
        DO J=1,3
           CALL SHAPE(XI(I),XI(J))
           CALL BMATRX(X,Y,BMAT,DET,R)
           CALL MULT1(DMAT,BMAT,BDB,C,4,8,4)

           DET = W(I)*W(J)*DET*R      

           S(1:8,1:8) = S(1:8,1:8) + BDB(1:8,1:8)*DET
           CN(1:7:2)  = CN(1:7:2) + (BMAT(1,1:7:2)+BMAT(4,1:7:2))*DET
           CN(2:8:2)  = CN(2:8:2) + (BMAT(2,2:8:2))*DET
           ECC(1:4)   = ECC(1:4)  +  AN(1:4)*DET	!  Lumped Mass Matrix
        ENDDO    !J
      ENDDO    !I

      RETURN
      END SUBROUTINE STIFF

C..........................................................................

      SUBROUTINE UPDATE(XG,YG,ICO,IX,TDISP,Ps,dPT)
      USE GCONTROL; USE MPROPERTY; USE NUM_INT; USE MPOINT
      IMPLICIT NONE
      DOUBLE PRECISION :: X(4),Y(4),XG(*),YG(*),U(8)
      DOUBLE PRECISION :: TDISP(*), Ps(*), Ps0, dPT(*),DS(4)
      DOUBLE PRECISION :: R,DET,DMAT(4,4),BMAT(4,8),STRAN(4)
      INTEGER :: IEL,J,K, IS,ico(5,1),IX(*),LJ(8)

C**   UPDATE STRESSES

      DO IEL=1,NEL
        CALL LOCAL(XG,YG,X,Y,IS,ICO,IEL,LJ,IX)
        CALL GLOLOC(SXX,SYY,SZZ,SXY,DS,IEL,1)

        U(1:8) = 0.D0
        DO J=1,8
          IF(LJ(J) > 0) U(J)=TDISP(LJ(J))
        ENDDO

        CALL DMATRX(EMOD(IS),ANV(IS),DMAT)
        CALL SHAPE(0.d0,0.d0)
        CALL BMATRX(X,Y,BMAT,DET,R) 

C** CALCULATE THE STRAIN INCREMENTS

        DO K=1,4
          STRAN(K)  = DOT_PRODUCT(BMAT(K,:),U(:))
        ENDDO

C**  CALCULATE THE STRESSES AT THE CENTROID

        DO K = 1,4
          DS(K) = DOT_PRODUCT(DMAT(K,:),STRAN(:))
        ENDDO

        CALL GLOLOC(SXX,SYY,SZZ,SXY,DS,IEL,0) ! effective stress saved
        Ps0 = Ps(IEL)                         ! initial Ps 
        Ps(IEL)  = -(SXX(IEL)+SYY(IEL)+SZZ(IEL))/3.d0 ! update  Ps
        dPT(IEL) = dPT(IEL) + (Ps(IEL)-Ps0)     ! update dPT	       
      ENDDO		! IEL							! 2nd stage update

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
      END	SUBROUTINE GLOLOC
	  
C...........................................................................

      SUBROUTINE LOCAL(XG,X,ICO,IEL,LJ,IX)
      IMPLICIT NONE
      DOUBLE PRECISION :: XG(*),X(2)
      INTEGER ::  IEL,LJ(*),IX(*),ico(2,*)

      X(1:2)=XG(ICO(1:2,IEL))
      LJ(2)=IX(ICO(1:2,IEL))

      RETURN
      END SUBROUTINE LOCAL
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

C..................................................................

      SUBROUTINE MatrixXa(B,X,S,NNET,LBAND) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(1), B(1), X(1)
      B(1:NNET) = 0.d0

C **  NOT USED FOR GLOBAL MULTIPLICATION IN THIS CODE
   
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

