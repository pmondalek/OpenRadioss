Chd|====================================================================
Chd|  SIGEPS57C                     source/materials/mat/mat057/sigeps57c.F
Chd|-- called by -----------
Chd|        MULAWC                        source/materials/mat_share/mulawc.F
Chd|-- calls ---------------
Chd|        VINTER                        source/tools/curve/vinter.F   
Chd|        FINTER                        source/tools/curve/finter.F   
Chd|====================================================================
      SUBROUTINE SIGEPS57C(
     1     NEL    ,NUPARAM,NUVAR   ,NFUNC   ,IFUNC  ,
     2     NPF    ,NPT    ,IPT     ,
     2     TF     ,TIME   ,TIMESTEP,UPARAM  ,RHO0   ,
     3     THKLY   ,
     4     EPSPXX ,EPSPYY ,EPSPXY  ,EPSPYZ  ,EPSPZX ,
     5     DEPSXX ,DEPSYY ,DEPSXY  ,DEPSYZ  ,DEPSZX ,
     6     EPSXX  ,EPSYY  ,EPSXY   ,EPSYZ   ,EPSZX  ,
     7     SIGOXX ,SIGOYY ,SIGOXY  ,SIGOYZ  ,SIGOZX ,
     8     SIGNXX ,SIGNYY ,SIGNXY  ,SIGNYZ  ,SIGNZX ,
     A     SOUNDSP,VISCMAX,THK     ,PLA     ,UVAR   ,
     B     OFF    ,NGL    ,SHF     ,ETSE    ,SIGY   ,
     C     HARDM  ,SEQ_OUTPUT,ISRATE,EPSP   ,INLOC  ,
     D     DPLANL ,SIGB)
Chd|====================================================================
C=======================================================================
C   BARLAT 89 material model
C=======================================================================
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "com01_c.inc"
#include      "sms_c.inc"
C---------+---------+---+---+--------------------------------------------
C VAR     | SIZE    |TYP| RW| DEFINITION
C---------+---------+---+---+--------------------------------------------
C NEL    |  1      | I | R | SIZE OF THE ELEMENT GROUP NEL F
C NUPARAM |  1      | I | R | SIZE OF THE USER PARAMETER ARRAY
C NUVAR   |  1      | I | R | NUMBER OF USER ELEMENT VARIABLES
C---------+---------+---+---+--------------------------------------------
C NFUNC   |  1      | I | R | NUMBER FUNCTION USED FOR THIS USER LAW
C IFUNC   | NFUNC   | I | R | FUNCTION INDEX 
C NPF     |  *      | I | R | FUNCTION ARRAY   
C IFLAG   |  *      | I | R | GEOMETRICAL FLAGS   
C TF      |  *      | F | R | FUNCTION ARRAY 
C---------+---------+---+---+--------------------------------------------
C TIME    |  1      | F | R | CURRENT TIME
C TIMESTEP|  1      | F | R | CURRENT TIME STEP
C UPARAM  | NUPARAM | F | R | USER MATERIAL PARAMETER ARRAY
C RHO0    | NEL    | F | R | INITIAL DENSITY
C THKLY   | NEL    | F | R | LAYER THICKNESS
C EPSPXX  | NEL    | F | R | STRAIN RATE XX
C EPSPYY  | NEL    | F | R | STRAIN RATE YY
C ...     |         |   |   |
C DEPSXX  | NEL    | F | R | STRAIN INCREMENT XX
C DEPSYY  | NEL    | F | R | STRAIN INCREMENT YY
C ...     |         |   |   |
C SIGOXX  | NEL    | F | R | OLD ELASTO PLASTIC STRESS XX 
C SIGOYY  | NEL    | F | R | OLD ELASTO PLASTIC STRESS YY
C ...     |         |   |   |    
C---------+---------+---+---+--------------------------------------------
C SIGNXX  | NEL    | F | W | NEW ELASTO PLASTIC STRESS XX
C SIGNYY  | NEL    | F | W | NEW ELASTO PLASTIC STRESS YY
C ...     |         |   |   |
C SOUNDSP | NEL    | F | W | SOUND SPEED (NEEDED FOR TIME STEP)
C---------+---------+---+---+--------------------------------------------
C THK     | NEL    | F |R/W| THICKNESS
C PLA     | NEL    | F |R/W| PLASTIC STRAIN
C UVAR    |NEL*NUVAR| F |R/W| USER ELEMENT VARIABLE ARRAY
C OFF     | NEL    | F |R/W| DELETED ELEMENT FLAG (=1. ON, =0. OFF)
C---------+---------+---+---+--------------------------------------------
C-----------------------------------------------
C   I N P U T   A r g u m e n t s
C-----------------------------------------------
      INTEGER,INTENT(IN)        :: NEL, NUVAR, NPT, IPT, NUPARAM,ISRATE,INLOC 
      INTEGER,DIMENSION(NEL),INTENT(IN) ::    NGL 
      my_real ,INTENT(IN) :: TIME,TIMESTEP
      my_real ,DIMENSION(NUPARAM),INTENT(IN) :: UPARAM
      my_real ,DIMENSION(NEL) ,INTENT(IN) :: EPSPXX ,EPSPYY ,EPSPXY ,EPSPYZ ,EPSPZX , 
     .   DEPSXX ,DEPSYY ,DEPSXY ,DEPSYZ ,DEPSZX , 
     .   EPSXX  ,EPSYY  , EPSXY  ,EPSYZ  ,EPSZX  ,
     .   SIGOXX ,SIGOYY ,SIGOXY ,SIGOYZ ,SIGOZX ,
     .   RHO0,SHF,THKLY,DPLANL 
      my_real ,DIMENSION(NEL)      ,INTENT(INOUT) :: PLA, OFF,THK
      my_real ,DIMENSION(NEL,3)    ,INTENT(INOUT) :: SIGB
      my_real ,DIMENSION(NEL,NUVAR),INTENT(INOUT) :: UVAR 
      my_real ,DIMENSION(NEL)      ,INTENT(OUT)   :: SIGNXX ,SIGNYY ,SIGNXY ,SIGNYZ ,SIGNZX , 
     .    SOUNDSP ,VISCMAX ,ETSE ,SIGY ,SEQ_OUTPUT ,HARDM 
      my_real ,DIMENSION(NEL)   :: EPSP 
C-----------------------------------------------
C   VARIABLES FOR FUNCTION INTERPOLATION 
C-----------------------------------------------
      INTEGER NPF(*), NFUNC, IFUNC(NFUNC)
      my_real FINTER ,TF(*)
      EXTERNAL FINTER
C-----------------------------------------------
C        Y = FINTER(IFUNC(J),X,NPF,TF,DYDX)
C        Y       : y = f(x)
C        X       : x
C        DYDX    : f'(x) = dy/dx
C        IFUNC(J): FUNCTION INDEX
C              J : FIRST(J=1), SECOND(J=2) .. FUNCTION USED FOR THIS LAW
C        NPF,TF  : FUNCTION PARAMETER
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER I,II,J,NRATE,J1,J2,N,NINDX,ITER,NITER,NS,
     .        NPARAM,IDFUNCE,NINC, OPTE,NMAX,IFLAGINT,INTM,
     .        IAD1(NEL),IPOS1(NEL),ILEN1(NEL),
     .        IAD2(NEL),IPOS2(NEL),ILEN2(NEL),
     .        JJ(NEL),INDX(NEL)

      my_real A, DEZZ, EPST,NU,ONE_M,A01, A02, A03,M, 
     .        A12,FISOKIN,NORM1,EPSMAX,EPSR1 ,EPSR2 ,CE,EINF 

      my_real SEFF,DK1_DSIGXX,DK1_DSIGYY,DK2_DSIGXX,DK2_DSIGYY,DF_DK1,NORM2,
     .        DF_DK2,NORMXX,NORMYY,DK2_DSIGXY,NORMXY,SIG_DFDSIG, 
     .        DDLAM,DYLD_DPLA,DF_DLAM,DFDSIG_N, DPLA_DLAM,DYLD_DLAM,
     .        DSXX,DSYY,DSXY,DEXX,DEYY,DEXY,ALPHA
      my_real , DIMENSION(NEL,2) :: YFAC, RATE
      my_real , DIMENSION(NEL)   ::E , A1 ,A2 ,G , G3 ,           
     .        DYDX1 ,DYDX2 ,SVM ,NORM,
     .        YLD  , Y1 , Y2 ,FRATE ,
     .        FAIL , H  ,GS  ,HK  ,        
     .        ESCALE ,C1 ,DYDXE ,YLDN ,SIGBXX,SIGBYY,SIGBXY,
     .        SIGEXX,SIGEYY,SIGEXY

      my_real, DIMENSION(NEL) :: K1,K2,DEPLZZ,DEELZZ,DPLA,DEPLXX,DEPLYY,DEPLXY,
     .                           YLD0,PHI,DPHI,DDEP,DYDX,DEPSZZ
C-----------------------------------------------
        NITER = 3
        NS = 10
        IFLAGINT = 2
        A01 = UPARAM(7)
        A02 = UPARAM(8)!!!A02 = C*2^M  !!!!
        A03 = UPARAM(9)
        A12 = UPARAM(10)
        NU  = UPARAM(6)           
        NRATE  = NINT(UPARAM(1))
        EPSMAX = UPARAM(2*NRATE+11)
        EPSR1  = UPARAM(2*NRATE+12)
        EPSR2  = UPARAM(2*NRATE+13)
        M      = UPARAM(2*NRATE+14)
        FISOKIN= UPARAM(2*NRATE+18)
        OPTE   = UPARAM(2*NRATE+20)
        EINF   = UPARAM(2*NRATE+21)
        DO I=1,NEL
          E(I)   = UPARAM(2)
          A1(I)  = UPARAM(3)
          A2(I)  = UPARAM(4)
          G(I)   = UPARAM(5)
          GS(I)  = G(I)*SHF(I)
        ENDDO
        CE     = UPARAM(2*NRATE+22)

        ONE_M =ONE/M
        IF (MOD(M,ONE)==ZERO)IFLAGINT = 1 !M is integer
        IF(IFLAGINT == 1) INTM = NINT(M)

C       
       IF(EPSMAX==ZERO)THEN
         IF(TF(NPF(IFUNC(1)+1)-1)==ZERO)THEN
          EPSMAX=TF(NPF(IFUNC(1)+1)-2)
         ELSE
          EPSMAX= EP30
         ENDIF
       ENDIF
       ! young modulus dependency on plastic strain
       IF (OPTE == 1)THEN   
         IDFUNCE = UPARAM(2*NRATE+19)
         DO I=1,NEL    
           IF(PLA(I) > ZERO)THEN                                                        
            ESCALE(I) = FINTER(IFUNC(IDFUNCE),PLA(I),NPF,TF,DYDXE(I))   
           ENDIF
         ENDDO  
         DO I=1,NEL    
           IF(PLA(I) > ZERO)THEN 
            E(I) =  ESCALE(I)* E(I)   
            A1(I) = E(I)/(ONE - NU *NU )
            A2(I) = NU *A1(I) 
            G(I) =  HALF*E(I)/(ONE+NU ) 
            GS(I)   =  G(I)*SHF(I)                                                 
           ENDIF  
         ENDDO                                                                                      
       ELSEIF ( CE /= ZERO) THEN      
         DO I=1,NEL                                   
          IF(PLA(I) > ZERO)THEN                                                        
            E(I) = E(I)-(E(I)-EINF)*(ONE-EXP(-CE*PLA(I))) 
            A1(I) = E(I)/(ONE - NU *NU )
            A2(I) = NU *A1(I)                                             
            G(I)  =  HALF*E(I)/(ONE+NU )   
            GS(I) =  G(I)*SHF(I)                                                                              
          ENDIF     
         ENDDO 
       ENDIF
c-----------------------------------------------
c     COMPUTATION OF THE TRIAL STRESS TENSOR
c-----------------------------------------------
       ! Computation of the trial stress tensor

       IF (FISOKIN > ZERO) THEN 
         DO I=1,NEL
          SIGNXX(I) = SIGOXX(I) - SIGB(I,1) + A1(I)*DEPSXX(I) + A2(I)*DEPSYY(I)
          SIGNYY(I) = SIGOYY(I) - SIGB(I,2) + A2(I)*DEPSXX(I) + A1(I)*DEPSYY(I)
          SIGNXY(I) = SIGOXY(I) - SIGB(I,3) + G(I)*DEPSXY(I)  
        ENDDO
      ELSE
        DO I=1,NEL
          SIGNXX(I)  = SIGOXX(I) + A1(I)*DEPSXX(I) + A2(I)*DEPSYY(I)
          SIGNYY(I)  = SIGOYY(I) + A2(I)*DEPSXX(I) + A1(I)*DEPSYY(I)
          SIGNXY(I)  = SIGOXY(I) + G(I) *DEPSXY(I)
        ENDDO
      ENDIF
      SIGEXX(1:NEL) = SIGNXX(1:NEL)
      SIGEYY(1:NEL) = SIGNYY(1:NEL)
      SIGEXY(1:NEL) = SIGNXY(1:NEL)
      DO I=1,NEL
        SIGNYZ(I)  = SIGOYZ(I) + GS(I)*DEPSYZ(I)
        SIGNZX(I)  = SIGOZX(I) + GS(I)*DEPSZX(I)
        ! Sound speed and tangent modulus
        SOUNDSP(I) = SQRT(A1(I)/RHO0(I))
        ETSE(I)    = ONE
        DPLA(I)    = ZERO
        DEPSZZ(I)  = ZERO
        DEPLZZ(I)  = ZERO
      ENDDO
      !normalisation 
      
      DO I=1,NEL
        ! Trial stress tensor
        NORM(I)   =  UPARAM(2)/EP03
        SIGNXX(I) =  SIGNXX(I)/NORM(I) 
        SIGNYY(I) =  SIGNYY(I)/NORM(I) 
        SIGNXY(I) =  SIGNXY(I)/NORM(I) 
      ENDDO
C-------------------
C     FAILURE STRAIN 
C-------------------
      DO I=1,NEL
        EPST   = HALF*( EPSXX(I)+EPSYY(I)
     .          + SQRT( (EPSXX(I)-EPSYY(I))*(EPSXX(I)-EPSYY(I))
     .                 + EPSXY(I)*EPSXY(I) ) )
        FAIL(I) = MAX(ZERO,MIN(ONE,(EPSR2-EPST)/(EPSR2-EPSR1)))        
      ENDDO                  
c-----------------------------------------------
c     COMPUTATION OF STRAIN-RATE
c-----------------------------------------------
      IF (ISRATE == 0) THEN
        DO I=1,NEL
          EPSP(I) = HALF*( ABS(EPSPXX(I)+EPSPYY(I))
     .                + SQRT( (EPSPXX(I)-EPSPYY(I))*(EPSPXX(I)-EPSPYY(I))
     .                + EPSPXY(I)*EPSPXY(I) ) )
        ENDDO
      ENDIF
C------------------------------------------------
          ! Computation of the yield stress
C------------------------------------------------
      DO I=1,NEL
         JJ(I) = 1
         DO J=2,NRATE-1
          IF( EPSP(I) >= UPARAM(NS+J) ) JJ(I) = J
         ENDDO
      ENDDO

      DO I=1,NEL
        RATE(I,1)=UPARAM(NS+JJ(I))
        RATE(I,2)=UPARAM(NS+JJ(I)+1)
        FRATE(I) = (EPSP(I) - RATE(I,1))/(RATE(I,2) - RATE(I,1))
        YFAC(I,1)=UPARAM(NS+NRATE+JJ(I))
        YFAC(I,2)=UPARAM(NS+NRATE+JJ(I)+1)
      ENDDO

      DO I=1,NEL
        J1 = JJ(I)
        J2 = J1+1
        IPOS1(I) = NINT(UVAR(I,J1+4))
        IAD1(I)  = NPF(IFUNC(J1)) / 2 + 1
        ILEN1(I) = NPF(IFUNC(J1)+1) / 2 - IAD1(I) - IPOS1(I)
        IPOS2(I) = NINT(UVAR(I,J2+4))
        IAD2(I)  = NPF(IFUNC(J2)) / 2 + 1
        ILEN2(I) = NPF(IFUNC(J2)+1) / 2 - IAD2(I) - IPOS2(I)
      ENDDO
C
      CALL VINTER(TF,IAD1,IPOS1,ILEN1,NEL,PLA,DYDX1,Y1)
      CALL VINTER(TF,IAD2,IPOS2,ILEN2,NEL,PLA,DYDX2,Y2)
C

      DO I=1,NEL
        J1 = JJ(I)
        J2 = J1+1 
        Y1(I)=Y1(I)*YFAC(I,1)
        Y2(I)=Y2(I)*YFAC(I,2)
        YLD(I) = FAIL(I)*(Y1(I)    + FRATE(I)*(Y2(I)-Y1(I)))
        YLD(I) = MAX(YLD(I),EM20)
        DYDX1(I) = DYDX1(I)*YFAC(I,1)
        DYDX2(I) = DYDX2(I)*YFAC(I,2)
        DYDX(I) =(ONE - FISOKIN)*  FAIL(I)*(DYDX1(I) + FRATE(I)*(DYDX2(I)-DYDX1(I)) )
C         ECROUISSAGE CINEMATIQUE
        Y1(I)=TF(NPF(IFUNC(J1))+1)*YFAC(I,1)
        Y2(I)=TF(NPF(IFUNC(J2))+1)*YFAC(I,2)
        YLD0(I) = FAIL(I)*(Y1(I)  + FRATE(I) * ( Y2(I)-Y1(I) ))
        YLD(I) = (ONE -FISOKIN) * YLD(I) + FISOKIN * YLD0(I)      
        UVAR(I,J1+4) = IPOS1(I)
        UVAR(I,J2+4) = IPOS2(I)
      ENDDO      
c------------------------------
        ! Computation of the Barlat equivalent stress      
c------------------------------ 
      IF (IFLAGINT==1) THEN 
        DO I=1,NEL
          K1(I)   = (SIGNXX(I) + A03*SIGNYY(I)) /TWO
          K2(I)   = SQRT(((SIGNXX(I) - A03*SIGNYY(I))/TWO)**2 + (A12*SIGNXY(I))**2) 
          SVM(I)  = A01*ABS(K1(I)+K2(I))**INTM + A01*ABS(K1(I)-K2(I))**INTM + A02*(K2(I))**INTM
          IF (SVM(I)>ZERO )THEN
          SVM(I)  = EXP(ONE_M *LOG (HALF * SVM(I) ) )   
          ELSE 
            SVM(I)  = ZERO
          ENDIF 
          PHI(I)    = SVM(I)*NORM(I) - YLD(I) 
        ENDDO
      ELSE  

        DO I=1,NEL
          K1(I)   = (SIGNXX(I) + A03*SIGNYY(I)) /TWO
          K2(I)   = SQRT(((SIGNXX(I) - A03*SIGNYY(I))/TWO)**2 + (A12*SIGNXY(I))**2) 
          SVM(I) = A01* EXP( M *LOG (MAX( ABS(K1(I)+K2(I)),EM20 ) )) 
     .           + A01* EXP( M *LOG (MAX( ABS(K1(I)-K2(I)),EM20 ) )) 
     .           + A02* EXP( M *LOG (MAX( K2(I),EM20 ) ))
          IF (SVM(I)>ZERO )THEN
            SVM(I)  = EXP(ONE_M *LOG (HALF * SVM(I) ) )   
          ELSE 
            SVM(I)  = ZERO
          ENDIF 
          PHI(I)    = SVM(I)*NORM(I) - YLD(I) 
        ENDDO      
      ENDIF 
c-----------------------------------------------
c       CHECKING THE YIELD CONDITION
c-----------------------------------------------
      NINDX  = 0
      DO I=1,NEL
        IF (PHI(I) > ZERO .AND.OFF(I) == ONE)THEN ! Plastic Loading
          NINDX = NINDX + 1
          INDX(NINDX)  = I 
        ENDIF
      ENDDO
c-----------------------------------------------
c         COMPUTING THE RETURN MAPPING
c-----------------------------------------------      
      IF (IFLAGINT==1) THEN 
#include "vectorize.inc" 
       DO II=1,NINDX                                              
        I = INDX(II) 
        J1 = JJ(I)
        J2 = J1+1
        DEPLXX(I) = ZERO
        DEPLYY(I) = ZERO
        DEPLXY(I) = ZERO      
        DO ITER = 1,NITER
            !  Computation of the derivatives
          DF_DK1 = (  EXP((ONE-M)*LOG  (SVM(I)))    ) *(A01/TWO)*(
     .                           +  SIGN(ONE,K1(I)+K2(I))*(ABS(K1(I)+K2(I))**(INTM-1)) 
     .                           +  SIGN(ONE,K1(I)-K2(I))*(ABS(K1(I)-K2(I))**(INTM-1)))
  
          DF_DK2 = (  EXP((ONE-M)*LOG  (SVM(I)))     ) *( (A01/TWO) * (
     .                           +  SIGN(ONE,K1(I)+K2(I))*(ABS(K1(I)+K2(I))**(INTM-1)) 
     .                           -  SIGN(ONE,K1(I)-K2(I))*(ABS(K1(I)-K2(I))**(INTM-1)) ) 
     .                           +  A02/TWO*((K2(I))**(INTM-1))   )
          DK1_DSIGXX = HALF 
          DK1_DSIGYY = A03/TWO
          DK2_DSIGXX = (SIGNXX(I)-A03*SIGNYY(I)) /(MAX(FOUR*K2(I),EM20))
          DK2_DSIGYY = -A03*(SIGNXX(I)-A03*SIGNYY(I))/(MAX(FOUR*K2(I),EM20))
          DK2_DSIGXY = (A12**2)*SIGNXY(I)/MAX(K2(I),EM20)
          !  Assembling the normal 
          NORMXX     = DF_DK1*DK1_DSIGXX+ DF_DK2*DK2_DSIGXX
          NORMYY     = DF_DK1*DK1_DSIGYY+ DF_DK2*DK2_DSIGYY
          NORMXY     = DF_DK2*DK2_DSIGXY


          ! df/ddlam = DF/DSIG * DSIG/Ddlam ( = df_dsig = N *Cel * N)
          DF_DLAM  = NORMXX * (A1(I)*NORMXX  +  A2(I)*NORMYY)
     .          +    NORMYY * (A1(I)*NORMYY  +  A2(I)*NORMXX)
     .          +    NORMXY * NORMXY * G(I)

          !compute derivative of effective plastic strain vs dlam
          DFDSIG_N   =  SIGNXX(I) * NORMXX
     .                + SIGNYY(I) * NORMYY
     .                + SIGNXY(I) * NORMXY 
          DPLA_DLAM  = DFDSIG_N *NORM(I)/ YLD(I) 
          !derivative of YLD to dlam
          DYLD_DPLA  = DYDX(I) 
          DYLD_DLAM  = DYLD_DPLA * DPLA_DLAM
          !d phi/ddlam
          DPHI(I) = - DF_DLAM - DYLD_DLAM
          DPHI(I) = SIGN(MAX(ABS(DPHI(I)),EM20) ,DPHI(I)) 

           ! Plastic multiplier increment
          DDLAM   = -PHI(I)/DPHI(I)     
          
          ! Plastic strains tensor update              
          DEPLXX(I) = DDLAM*NORMXX
          DEPLYY(I) = DDLAM*NORMYY
          DEPLXY(I) = DDLAM*NORMXY

c
          ! Cumulated plastic strain update
          DDEP(I)   = DDLAM * DPLA_DLAM
          PLA(I)    = PLA(I) + DDEP(I)

c
          ! Total plastic strain increment on the time step
          DPLA(I) = DPLA(I) + DDEP(I) 
          DPLA(I) = MAX(DPLA(I),ZERO)
c
          ! Cauchy stress tensor update              
          SIGNXX(I) = SIGNXX(I)*NORM(I) - A1(I)*DEPLXX(I)  - A2(I)*DEPLYY(I)                
          SIGNYY(I) = SIGNYY(I)*NORM(I) - A2(I)*DEPLXX(I)  - A1(I)*DEPLYY(I)            
          SIGNXY(I) = SIGNXY(I)*NORM(I) - G(I)*DEPLXY(I)  

          SIGNXX(I) = SIGNXX(I)/NORM(I) 
          SIGNYY(I) = SIGNYY(I)/NORM(I) 
          SIGNXY(I) = SIGNXY(I)/NORM(I) 

          K1(I)      = (SIGNXX(I) + A03*SIGNYY(I))/TWO
          K2(I)      = SQRT(((SIGNXX(I) - A03*SIGNYY(I))/TWO)**2 + (A12*SIGNXY(I))**2)
          SVM(I)     = A01*ABS(K1(I)+K2(I))**INTM + A01*ABS(K1(I)-K2(I))**INTM+A02*(K2(I))**INTM
          SVM(I)     = EXP(ONE_M *LOG(HALF*MAX(SVM(I) ,EM20))) 

          Y1(I)  = YFAC(I,1)*FINTER(IFUNC(J1),PLA(I),NPF,TF,DYDX1(I)) 
          Y2(I)  = YFAC(I,2)*FINTER(IFUNC(J2),PLA(I),NPF,TF,DYDX2(I)) 
          DYDX1(I) = DYDX1(I)*YFAC(I,1)
          DYDX2(I) = DYDX2(I)*YFAC(I,2)
  
          YLD(I) = (Y1(I)    + FRATE(I)*(Y2(I)-Y1(I)))
          YLD(I) = MAX(YLD(I),EM20)
          YLD(I) = (ONE - FISOKIN) * YLD(I) + FISOKIN * YLD0(I)
          DYDX(I)= (ONE - FISOKIN)* FAIL(I)*(DYDX1(I) + FRATE(I)*(DYDX2(I)-DYDX1(I)))

          
          PHI(I) = SVM(I)*NORM(I) - YLD(I)  

          ! Thickness plastic strain update
          DEPLZZ(I) = DEPLZZ(I) - (DEPLXX(I)+DEPLYY(I))

        ENDDO
       ENDDO
      ELSE ! M is real value
#include "vectorize.inc" 
        DO II=1,NINDX                                              
         I = INDX(II) 
         J1 = JJ(I)
         J2 = J1+1
         DEPLXX(I) = ZERO
         DEPLYY(I) = ZERO
         DEPLXY(I) = ZERO      
         DO ITER = 1,NITER
             !  Computation of the derivatives
           DF_DK1 = (  EXP((ONE-M)*LOG  (SVM(I)))    ) *(A01/TWO) * (
     .         SIGN(ONE,K1(I)+K2(I))* ( EXP( (M-ONE) * LOG ( MAX(ABS(K1(I)+K2(I)) , EM20)  ) ) )
     .      +  SIGN(ONE,K1(I)-K2(I))* ( EXP( (M-ONE) * LOG ( MAX(ABS(K1(I)-K2(I)) , EM20)  ) ) ) )
   
           DF_DK2 = (  EXP((ONE-M)*LOG  (SVM(I)))     ) *( (A01/TWO)*(
     .          SIGN(ONE,K1(I)+K2(I))* (EXP((M-ONE) * LOG ( MAX(ABS(K1(I)+K2(I)) , EM20)  ) ) )
     .       -  SIGN(ONE,K1(I)-K2(I))* (EXP((M-ONE) * LOG ( MAX(ABS(K1(I)-K2(I)) , EM20)  ) ) ) )
     .       +  A02/TWO*               (EXP((M-ONE) * LOG ( MAX(K2(I), EM20)   ) )) )
           DK1_DSIGXX = HALF 
           DK1_DSIGYY = A03/TWO
           DK2_DSIGXX = (SIGNXX(I)-A03*SIGNYY(I)) /(MAX(FOUR*K2(I),EM20))
           DK2_DSIGYY = -A03*(SIGNXX(I)-A03*SIGNYY(I))/(MAX(FOUR*K2(I),EM20))
           DK2_DSIGXY = (A12**2)*SIGNXY(I)/MAX(K2(I),EM20)
           !  Assembling the normal 
           NORMXX     = DF_DK1*DK1_DSIGXX+ DF_DK2*DK2_DSIGXX
           NORMYY     = DF_DK1*DK1_DSIGYY+ DF_DK2*DK2_DSIGYY
           NORMXY     = DF_DK2*DK2_DSIGXY
 
 
           ! df/ddlam = DF/DSIG * DSIG/Ddlam ( = df_dsig = N *Cel * N)
           DF_DLAM  = NORMXX * (A1(I)*NORMXX  +  A2(I)*NORMYY)
     .          +    NORMYY * (A1(I)*NORMYY  +  A2(I)*NORMXX)
     .          +    NORMXY * NORMXY * G(I)
 
           !compute derivative of effective plastic strain vs dlam
           DFDSIG_N   =  SIGNXX(I) * NORMXX
     .                + SIGNYY(I) * NORMYY
     .                + SIGNXY(I) * NORMXY 
           DPLA_DLAM  = DFDSIG_N *NORM(I)/ YLD(I) 
           !derivative of YLD to dlam
           DYLD_DPLA  = DYDX(I) 
           DYLD_DLAM  = DYLD_DPLA * DPLA_DLAM
           !d phi/ddlam
           DPHI(I) = - DF_DLAM - DYLD_DLAM
           DPHI(I) = SIGN(MAX(ABS(DPHI(I)),EM20) ,DPHI(I)) 
 
            ! Plastic multiplier increment
           DDLAM   = -PHI(I)/DPHI(I)     
           
           ! Plastic strains tensor update              
           DEPLXX(I) = DDLAM*NORMXX
           DEPLYY(I) = DDLAM*NORMYY
           DEPLXY(I) = DDLAM*NORMXY
 
c
           ! Cumulated plastic strain update
           DDEP(I)   = DDLAM * DPLA_DLAM
           PLA(I)    = PLA(I) + DDEP(I)
 
c
           ! Total plastic strain increment on the time step
           DPLA(I) = DPLA(I) + DDEP(I) 
           DPLA(I) = MAX(DPLA(I),ZERO)
c
           ! Cauchy stress tensor update              
           SIGNXX(I) = SIGNXX(I)*NORM(I) - A1(I)*DEPLXX(I)  - A2(I)*DEPLYY(I)                
           SIGNYY(I) = SIGNYY(I)*NORM(I) - A2(I)*DEPLXX(I)  - A1(I)*DEPLYY(I)            
           SIGNXY(I) = SIGNXY(I)*NORM(I) - G(I)*DEPLXY(I)  
 
           SIGNXX(I) = SIGNXX(I)/NORM(I) 
           SIGNYY(I) = SIGNYY(I)/NORM(I) 
           SIGNXY(I) = SIGNXY(I)/NORM(I) 
 
           K1(I)      = (SIGNXX(I) + A03*SIGNYY(I))/TWO
           K2(I)      = SQRT(((SIGNXX(I) - A03*SIGNYY(I))/TWO)**2 + (A12*SIGNXY(I))**2)
           SVM(I)     =  A01* EXP( M *LOG (MAX( ABS(K1(I)+K2(I)),EM20 ) )) 
     .                 + A01* EXP( M *LOG (MAX( ABS(K1(I)-K2(I)),EM20 ) )) 
     .                 + A02* EXP( M *LOG (MAX( K2(I),EM20 ) ))
                 
           SVM(I)     = EXP(ONE_M *LOG(HALF*MAX(SVM(I) ,EM20))) 
 
           Y1(I)  = YFAC(I,1)*FINTER(IFUNC(J1),PLA(I),NPF,TF,DYDX1(I)) 
           Y2(I)  = YFAC(I,2)*FINTER(IFUNC(J2),PLA(I),NPF,TF,DYDX2(I)) 
           DYDX1(I) = DYDX1(I)*YFAC(I,1)
           DYDX2(I) = DYDX2(I)*YFAC(I,2)
   
           YLD(I) = (Y1(I)    + FRATE(I)*(Y2(I)-Y1(I)))
           YLD(I) = MAX(YLD(I),EM20)
           YLD(I) = (ONE - FISOKIN) * YLD(I) + FISOKIN * YLD0(I)
           DYDX(I)= (ONE - FISOKIN)* FAIL(I)*(DYDX1(I) + FRATE(I)*(DYDX2(I)-DYDX1(I)))
 
           
           PHI(I) = SVM(I)*NORM(I) - YLD(I)  
 
           ! Thickness plastic strain update
           DEPLZZ(I) = DEPLZZ(I) - (DEPLXX(I)+DEPLYY(I))
 
         ENDDO
        ENDDO        
      ENDIF
c
      DO I=1,NEL
        IF(PLA(I)>EPSMAX.AND.OFF(I)==ONE)OFF(I)=ZEP80
      ENDDO
      !remove normalisation 
      DO I=1,NEL
        SIGNXX(I) = SIGNXX(I)*NORM(I) 
        SIGNYY(I) = SIGNYY(I)*NORM(I) 
        SIGNXY(I) = SIGNXY(I)*NORM(I) 
        SVM(I) =  SVM(I)*NORM(I)
 
      ENDDO      
      DO I=1,NEL
        IF(DPLA(I) > ZERO.AND.DYDX(I) >ZERO) THEN
          H(I)    = DYDX(I)
          ETSE(I) = H(I)/(H(I)+E(I))
        ELSE
          ETSE(I) = ONE 
        ENDIF             
        SEQ_OUTPUT(I) = SVM(I) 
        DEELZZ(I ) = -NU*(SIGNXX(I)-SIGOXX(I)+SIGNYY(I)-SIGOYY(I))/E(I)
      ENDDO
c      update bakstress      
      IF (FISOKIN > ZERO ) THEN
        DO I=1,NEL
C
          DSXX  = SIGEXX(I) - SIGNXX(I)                      
          DSYY  = SIGEYY(I) - SIGNYY(I)                      
          DSXY  = SIGEXY(I) - SIGNXY(I)                      
          DEXX  = (DSXX - NU*DSYY)                          
          DEYY  = (DSYY - NU*DSXX)                          
          DEXY  = TWO*(ONE+NU)*DSXY                         
          ALPHA = FISOKIN*H(I)/(E(I)+H(I)) * THIRD   
          SIGB(I,1) =  SIGB(I,1) + ALPHA*(FOUR*DEXX+TWO*DEYY)
          SIGB(I,2) =  SIGB(I,2) + ALPHA*(FOUR*DEYY+TWO*DEXX)
          SIGB(I,3) =  SIGB(I,3) + ALPHA*DEXY 
          SIGNXX(I) = SIGNXX(I) + SIGB(I,1)  
          SIGNYY(I) = SIGNYY(I) + SIGB(I,2)  
          SIGNXY(I) = SIGNXY(I) + SIGB(I,3)  
          UVAR(I,2) = SIGB(I,1)
          UVAR(I,3) = SIGB(I,2)
          UVAR(I,4) = SIGB(I,3)
       ENDDO
      ENDIF

!-------------------------------------------------------------
!     NON-LOCAL THICKNESS VARIATION
!-------------------------------------------------------------
      IF (INLOC > 0) THEN
        DO I = 1,NEL
          ! Computation of the normal to the yield surface
          K1(I)      = (SIGNXX(I) + A03*SIGNYY(I))/TWO
          K2(I)      = SQRT(((SIGNXX(I) - A03*SIGNYY(I))/TWO)**2 + (A12*SIGNXY(I))**2)
          SEFF    = A01*ABS(K1(I)+K2(I))**M + A01*ABS(K1(I)-K2(I))**M + A02*ABS(K2(I))**M
          SEFF    = EXP(ONE_M *LOG(HALF*MAX(SEFF ,EM20))) 
          
          DF_DK1  = (SEFF**(1-M))*(A01/TWO)*(
     .              +  SIGN(ONE,K1(I)+K2(I))*(ABS(K1(I)+K2(I))**(M-1)) 
     .              +  SIGN(ONE,K1(I)-K2(I))*(ABS(K1(I)-K2(I))**(M-1)))
          DF_DK2  = (SEFF**(1-M))*((A01/TWO)*(
     .              +  SIGN(ONE,K1(I)+K2(I))*(ABS(K1(I)+K2(I))**(M-1)) 
     .              -  SIGN(ONE,K1(I)-K2(I))*(ABS(K1(I)-K2(I))**(M-1))) 
     .              +  A02/TWO*(ABS(K2(I))**(M-1)))
          DK1_DSIGXX = HALF
          DK1_DSIGYY = A03/TWO
          DK2_DSIGXX = (SIGNXX(I)-A03*SIGNYY(I))/ (MAX(FOUR*K2(I),EM20))
          DK2_DSIGYY = -A03*(SIGNXX(I)-A03*SIGNYY(I))/(MAX(FOUR*K2(I),EM20))
          DK2_DSIGXY = (A12**2)*SIGNXY(I)/MAX(K2(I),EM20)
          NORMXX     = DF_DK1*DK1_DSIGXX + DF_DK2*DK2_DSIGXX
          NORMYY     = DF_DK1*DK1_DSIGYY + DF_DK2*DK2_DSIGYY
          NORMXY     = DF_DK2*DK2_DSIGXY
          SIG_DFDSIG = NORMXX*SIGNXX(I) + NORMYY*SIGNYY(I) + NORMXY*SIGNXY(I)
          ! Updating the thickness
          DEZZ       = -NU *((SIGNXX(I)-SIGOXX(I)+SIGNYY(I)-SIGOYY(I))/E(I))
          IF (SIG_DFDSIG /= ZERO) THEN
            DEZZ     = DEZZ - MAX(DPLANL(I),ZERO)*(SEFF/SIG_DFDSIG)*(NORMXX + NORMYY)
          ENDIF
          THK(I)     = THK(I) + DEZZ*THKLY(I)*OFF(I)   
        ENDDO  
      ELSE 
        DO I = 1,NEL
          DEPSZZ(I)  = DEELZZ(I) + DEPLZZ(I)
          THK(I)     = THK(I) + DEPSZZ(I)*THKLY(I)*OFF(I)
        ENDDO
      ENDIF
      HARDM(1:NEL) = DYDX(1:NEL)
      SIGY(1:NEL)  = YLD(1:NEL)
      !------------------------- 
      RETURN
      END  
