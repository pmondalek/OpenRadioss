Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2023 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
Chd|====================================================================
Chd|  SIGEPS190                      source/materials/mat/mat190/sigeps190.F
Chd|-- called by -----------
Chd|        MULAW                         source/materials/mat_share/mulaw.F
Chd|-- calls ---------------
Chd|        FINTER                        source/tools/curve/finter.F   
Chd|====================================================================
      SUBROUTINE SIGEPS190(
     1     NEL    ,NUPARAM,NUVAR   ,NFUNC   ,IFUNC   ,NPF    ,
     2     TF     ,TIME   ,TIMESTEP,UPARAM  ,RHO0    ,RHO    ,
     3     VOLUME ,EINT   ,
     4     EPSPXX ,EPSPYY ,EPSPZZ  ,EPSPXY  ,EPSPYZ  ,EPSPZX ,
     5     DEPSXX ,DEPSYY ,DEPSZZ  ,DEPSXY  ,DEPSYZ  ,DEPSZX ,
     6     EPSXX  ,EPSYY  ,EPSZZ   ,EPSXY   ,EPSYZ   ,EPSZX  ,
     7     SIGOXX ,SIGOYY ,SIGOZZ  ,SIGOXY  ,SIGOYZ  ,SIGOZX ,
     8     SIGNXX ,SIGNYY ,SIGNZZ  ,SIGNXY  ,SIGNYZ  ,SIGNZX ,
     9      MFXX   ,MFXY    ,MFXZ,MFYX, MFYY  , MFYZ  ,  
     a      MFZX   ,MFZY    ,MFZZ     ,             
     A      SOUNDSP,VISCMAX,UVAR    ,OFF     ,NGL     ,ISMSTR,
     B      EPSD  ,MATPARAM)
C User Material law for isotropic path dependent recoverable foam
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE TABLE_MOD
      USE INTERFACE_TABLE_MOD
      USE MATPARAM_DEF_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
#include      "comlock.inc"
#include      "scr05_c.inc"
#include      "impl1_c.inc"
#include      "tabsiz_c.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
C---------+---------+---+---+--------------------------------------------
C VAR     | SIZE    |TYP| RW| DEFINITION
C---------+---------+---+---+--------------------------------------------
C NEL     |  1      | I | R | SIZE OF THE ELEMENT GROUP NEL 
C NUPARAM |  1      | I | R | SIZE OF THE USER PARAMETER ARRAY
C NUVAR   |  1      | I | R | NUMBER OF USER ELEMENT VARIABLES
C---------+---------+---+---+--------------------------------------------
C MFUNC   |  1      | I | R | NUMBER FUNCTION USED FOR THIS USER LAW not used
C KFUNC   | NFUNC   | I | R | FUNCTION INDEX not used
C NPF     |  *      | I | R | FUNCTION ARRAY   
C TF      |  *      | F | R | FUNCTION ARRAY 
C---------+---------+---+---+--------------------------------------------
C TIME    |  1      | F | R | CURRENT TIME
C TIMESTEP|  1      | F | R | CURRENT TIME STEP
C UPARAM  | NUPARAM | F | R | USER MATERIAL PARAMETER ARRAY
C RHO0    | NEL     | F | R | INITIAL DENSITY
C RHO     | NEL     | F | R | DENSITY
C VOLUME  | NEL     | F | R | VOLUME
C EINT    | NEL     | F | R | TOTAL INTERNAL ENERGY
C EPSPXX  | NEL     | F | R | STRAIN RATE XX
C EPSPYY  | NEL     | F | R | STRAIN RATE YY
C ...     |         |   |   |
C DEPSXX  | NEL     | F | R | STRAIN INCREMENT XX
C DEPSYY  | NEL     | F | R | STRAIN INCREMENT YY
C ...     |         |   |   |
C EPSXX   | NEL     | F | R | STRAIN XX
C EPSYY   | NEL     | F | R | STRAIN YY
C ...     |         |   |   |
C SIG0XX  | NEL     | F | R | OLD ELASTO PLASTIC STRESS XX 
C SIG0YY  | NEL     | F | R | OLD ELASTO PLASTIC STRESS YY
C ...     |         |   |   |    
C---------+---------+---+---+--------------------------------------------
C SIGNXX  | NEL     | F | W | NEW ELASTO PLASTIC STRESS XX
C SIGNYY  | NEL     | F | W | NEW ELASTO PLASTIC STRESS YY
C ...     |         |   |   |
C SIGVXX  | NEL     | F | W | VISCOUS STRESS XX
C SIGVYY  | NEL     | F | W | VISCOUS STRESS YY
C ...     |         |   |   |
C SOUNDSP | NEL     | F | W | SOUND SPEED (NEEDED FOR TIME STEP)
C VISCMAX | NEL     | F | W | MAXIMUM DAMPING MODULUS(NEEDED FOR TIME STEP)
C---------+---------+---+---+--------------------------------------------
C UVAR    |NEL*NUVAR| F |R/W| USER ELEMENT VARIABLE ARRAY
C OFF     | NEL     | F |R/W| DELETED ELEMENT FLAG (=1. ON, =0. OFF)
C---------+---------+---+---+--------------------------------------------
C-------------------------------------------------------------------------
C-----------------------------------------------
C   I N P U T   A r g u m e n t s
C-----------------------------------------------
      INTEGER ,INTENT(IN) :: NEL,NFUNC,NUPARAM,NUVAR,ISMSTR
      INTEGER, DIMENSION(NEL),INTENT(IN) :: NGL
      my_real ,INTENT(IN) :: TIME,TIMESTEP
      INTEGER ,DIMENSION(NFUNC),INTENT(IN) :: IFUNC
      INTEGER ,DIMENSION(SNPC) ,INTENT(IN) :: NPF
      my_real ,DIMENSION(STF)  ,INTENT(IN) :: TF
      my_real ,DIMENSION(NUPARAM) ,INTENT(IN) :: UPARAM
      my_real ,DIMENSION(NEL) ,INTENT(IN) :: RHO,RHO0,VOLUME ,EINT,
     .         EPSXX,EPSYY,EPSZZ,EPSXY,EPSYZ,EPSZX,
     .         EPSPXX,EPSPYY,EPSPZZ,EPSPXY,EPSPYZ,EPSPZX,
     .         DEPSXX,DEPSYY,DEPSZZ,DEPSXY,DEPSYZ,DEPSZX,
     .         SIGOXX,SIGOYY,SIGOZZ,SIGOXY,SIGOYZ,SIGOZX,
     .         MFXX   ,MFXY   ,MFXZ ,MFYX, MFYY  , MFYZ  ,  
     .         MFZX   ,MFZY   ,MFZZ ,OFF                  
C-----------------------------------------------
C   O U T P U T   A r g u m e n t s
C-----------------------------------------------
      my_real ,DIMENSION(NEL) ,INTENT(OUT) ::
     .      SIGNXX ,SIGNYY ,SIGNZZ  ,SIGNXY  ,SIGNYZ  ,SIGNZX,VISCMAX
C-----------------------------------------------
C   I N P U T   O U T P U T   A r g u m e n t s 
C-----------------------------------------------
      my_real, DIMENSION(NEL) ,INTENT(INOUT) :: EPSD,SOUNDSP 
      my_real ,DIMENSION(NEL,NUVAR) ,INTENT(INOUT) :: UVAR
      TYPE(MATPARAM_STRUCT_) , TARGET :: MATPARAM
C----------------------------------------------------------------
c    VARIABLES FOR FUNCTION INTERPOLATION 
C----------------------------------------------------------------
      EXTERNAL FINTER
C----------------------------------------------------------------
C  L O C A L  V A R I B L E S
C----------------------------------------------------------------
      INTEGER I,IK,IJ 
      my_real, DIMENSION(NEL) :: 
        ! Deformation GRADIENT
     .   FPSXX,FPSYY,FPSZZ,
     .   FPSXY,FPSYX,FPSXZ,
     .   FPSZX,FPSYZ,FPSZY,
         ! Old Deformation GRADIENT
     .   FPOSXX,FPOSYY,FPOSZZ,
     .   FPOSXY,FPOSYX,FPOSXZ,
     .   FPOSZX,FPOSYZ,FPOSZY
      my_real
     .   G,K,NU,HU, SHAPE ,GS,KS,F(NEL,3,3),FO(NEL,3,3),EP1,
     .   EP2,EP3,EP4,EP5,EP6,ERT11,ERT12,ERT13,ERT21,
     .   ERT22,ERT23,ERT31,ERT32,ERT33,EPSP(3),strain(NEL,3),
     .   JAC(NEL),Jeq(NEL),B11(NEL),B22(NEL),B33(NEL),B12(NEL),
     .   B23(NEL),B13(NEL),Bsqr11(NEL),Bsqr22(NEL),I1bar(NEL),I2bar(NEL),
     .   Bsqr33(NEL),Bsqr12(NEL),Bsqr23(NEL),Bsqr13(NEL),H1(NEL),H2(NEL),H3(NEL),
     .   E,EMOD,
     .   ZXX(NEL),ZXY(NEL),ZYZ(NEL),ZXZ(NEL), ZYX(NEL), ZZY(NEL),
     .   ZZX(NEL),ZYY(NEL),ZZZ(NEL),WHYSMAX(NEL),
     .   ZXXOLD(NEL),ZXYOLD(NEL),ZYZOLD(NEL),ZXZOLD(NEL), 
     .   ZYXOLD(NEL), ZZYOLD(NEL),
     .   ZZXOLD(NEL),ZYYOLD(NEL),ZZZOLD(NEL),
     .   ZXXMID(NEL),ZXYMID(NEL),ZYZMID(NEL),ZXZMID(NEL), 
     .   ZYXMID(NEL), ZZYMID(NEL),
     .   ZZXMID(NEL),ZYYMID(NEL),ZZZMID(NEL),
     .   DZXX(NEL),DZXY(NEL),DZYZ(NEL),DZXZ(NEL), 
     .   DZYX(NEL), DZZY(NEL),
     .   DZZX(NEL),DZYY(NEL),DZZZ(NEL),
     .   spknodam(6,NEL),spknorate(6,NEL),spk(6,NEL),
     .   S1(NEL),S2(NEL),S3(NEL),S4(NEL),S5(NEL),S6(NEL),
     .   FS11(NEL),FS12(NEL),FS13(NEL),FS21(NEL),FS22(NEL),FS23(NEL),
     .   FS31(NEL),FS32(NEL),FS33(NEL),SCAL,SCALINT,SLOPEMAX(NEL),
     .   damage(NEL),DRATE(NEL,3),dfirst(3,6,NEL),cijkl(NEL,6,6),
     .   SLOPE(NEL,3),Zxxd(NEL),Zyyd(NEL),Zzzd(NEL),Zxyd(NEL),Zyzd(NEL),Zzxd(NEL),
     .   SIG(NEL,6), VAL(NEL,3), VEC(NEL,9),AV(6,NEL),DIRPRV(3,3,NEL),EV(NEL,3),EVV(3,NEL)
      INTEGER NUMTABL
      TYPE(TTABLE), DIMENSION(:) , POINTER :: TABLE
c===================================================================================
c      Isotropic path dependent recoverable foam law
c-------------------------------------------------------------------------
C      


       !print* ,'LOI 190 --------------------------------'
       EMOD = UPARAM(1)     
       E    = UPARAM(2)     
       K    = UPARAM(3)     
       NU   = UPARAM(4)     
       G    = UPARAM(5)     
       SCAL = UPARAM(6)     
       SCALINT = UPARAM(7)  
       HU      = UPARAM(8)       
       SHAPE   = UPARAM(9)    
       NUMTABL =  MATPARAM%NTABLE
       TABLE   => MATPARAM%TABLE(1:NUMTABL)
c===================================================================================
       !print*, 'NUMTABL ', NUMTABL
c     FPxx = Deformation Gradient
       DO I=1,NEL                  
         FPSXX(I) = ONE+MFXX(I)     
         FPSYY(I) = ONE+MFYY(I)     
         FPSZZ(I) = ONE+MFZZ(I)     
       END DO                      
       FPSXY(1:NEL) = MFXY(1:NEL)       
       FPSYX(1:NEL) = MFYX(1:NEL)       
       FPSXZ(1:NEL) = MFXZ(1:NEL)       
       FPSZX(1:NEL) = MFZX(1:NEL)       
       FPSYZ(1:NEL) = MFYZ(1:NEL)       
       FPSZY(1:NEL) = MFZY(1:NEL)       
C
c-------------------------------------------
c    FPxx = Old Deformation Gradient from previous timestep
c
       IF ( TIME == ZERO ) THEN             
         DO I=1,NEL                        
c         if (i.eq.1) print *,i,' uvar-0'  
           UVAR(I,1) = ONE	               
           UVAR(I,2) = ONE                
           UVAR(I,3) = ONE                 
           spknodam(1,i) = ZERO             
           spknodam(2,i) = ZERO             
           spknodam(3,i) = ZERO             
           spknodam(4,i) = ZERO             
           spknodam(5,i) = ZERO             
           spknodam(6,i) = ZERO             
           spknorate(1,i) = ZERO            
           spknorate(2,i) = ZERO            
           spknorate(3,i) = ZERO            
           spknorate(4,i) = ZERO            
           spknorate(5,i) = ZERO            
           spknorate(6,i) = ZERO            
         ENDDO                             
       ENDIF                               
      
c                              
       DO I=1,NEL                          
         FPOSXX(I) = UVAR(I,1)		           
         FPOSYY(I) = UVAR(I,2)             
         FPOSZZ(I) = UVAR(I,3)             
         FPOSXY(I) = UVAR(I,4)             
         FPOSYX(I) = UVAR(I,5)             
         FPOSXZ(I) = UVAR(I,6)             
         FPOSZX(I) = UVAR(I,7)             
         FPOSYZ(I) = UVAR(I,8)             
         FPOSZY(I) = UVAR(I,9)             
c
         UVAR(I,1) = FPSXX(I)              
         UVAR(I,2) = FPSYY(I)              
         UVAR(I,3) = FPSZZ(I)              
         UVAR(I,4) = FPSXY(I)              
         UVAR(I,5) = FPSYX(I)              
         UVAR(I,6) = FPSXZ(I)              
         UVAR(I,7) = FPSZX(I)              
         UVAR(I,8) = FPSYZ(I)              
         UVAR(I,9) = FPSZY(I)      
         spknorate(1,i) = UVAR(I,11)		          
         spknorate(2,i) = UVAR(I,12)            
         spknorate(3,i) = UVAR(I,13)            
         spknorate(4,i) = UVAR(I,14)            
         spknorate(5,i) = UVAR(I,15)            
         spknorate(6,i) = UVAR(I,16)            
         WHYSMAX(I)     = UVAR(I,17)        
       ENDDO                               
                             
c Path dependent solution
c-------------------------------------------
c  Loop over all elements
c
C         

       DO I=1,NEL
C
c      deformation gradient at t(n+1)
C
         F(I,1,1) = FPSXX(I)
         F(I,2,2) = FPSYY(I)
         F(I,3,3) = FPSZZ(I)
         F(I,1,2) = FPSXY(I)
         F(I,2,3) = FPSYZ(I)
         F(I,1,3) = FPSXZ(I)
         F(I,2,1) = FPSYX(I)
         F(I,3,2) = FPSZY(I)
         F(I,3,1) = FPSZX(I)
c
c        Green-Lagrange strains at Zt(n+1) = (FT * F - I)/2
c
         ZXX(i) = F(I,1,1)*F(I,1,1) + F(I,2,1)*F(I,2,1) + F(I,3,1)*F(I,3,1) 
         ZYY(i) = F(I,1,2)*F(I,1,2) + F(I,2,2)*F(I,2,2) + F(I,3,2)*F(I,3,2) 
         ZZZ(i) = F(I,1,3)*F(I,1,3) + F(I,2,3)*F(I,2,3) + F(I,3,3)*F(I,3,3) 
         ZXY(i) = F(I,1,1)*F(I,1,2) + F(I,2,1)*F(I,2,2) + F(I,3,1)*F(I,3,2) 
         ZYZ(i) = F(I,1,2)*F(I,1,3) + F(I,2,2)*F(I,2,3) + F(I,3,2)*F(I,3,3) 
         ZZX(i) = F(I,1,1)*F(I,1,3) + F(I,2,1)*F(I,2,3) + F(I,3,1)*F(I,3,3) 									
         ZXX(i) = (ZXX(i) - ONE) /TWO
         ZYY(i) = (ZYY(i) - ONE) /TWO
         ZZZ(i) = (ZZZ(i) - ONE) /TWO
         ZXY(i) = (ZXY(i)) /TWO
         ZYZ(i) = (ZYZ(i)) /TWO
         ZZX(i) = (ZZX(i)) /TWO
C
c        deformation gradient at t(n)
C
         FO(I,1,1) = FPOSXX(I)
         FO(I,2,2) = FPOSYY(I)
         FO(I,3,3) = FPOSZZ(I)
         FO(I,1,2) = FPOSXY(I)
         FO(I,2,3) = FPOSYZ(I)
         FO(I,1,3) = FPOSXZ(I)
         FO(I,2,1) = FPOSYX(I)
         FO(I,3,2) = FPOSZY(I)
         FO(I,3,1) = FPOSZX(I)
C
c        Green-Lagrange strains at t(n)
C
         ZXXOLD(i) = FO(I,1,1)*FO(I,1,1) + FO(I,2,1)*FO(I,2,1) + FO(I,3,1)*FO(I,3,1) 
         ZYYOLD(i) = FO(I,1,2)*FO(I,1,2) + FO(I,2,2)*FO(I,2,2) + FO(I,3,2)*FO(I,3,2) 
         ZZZOLD(i) = FO(I,1,3)*FO(I,1,3) + FO(I,2,3)*FO(I,2,3) + FO(I,3,3)*FO(I,3,3) 
         ZXYOLD(i) = FO(I,1,1)*FO(I,1,2) + FO(I,2,1)*FO(I,2,2) + FO(I,3,1)*FO(I,3,2) 
         ZYZOLD(i) = FO(I,1,2)*FO(I,1,3) + FO(I,2,2)*FO(I,2,3) + FO(I,3,2)*FO(I,3,3) 
         ZZXOLD(i) = FO(I,1,1)*FO(I,1,3) + FO(I,2,1)*FO(I,2,3) + FO(I,3,1)*FO(I,3,3) 									
         ZXXOLD(i)=(ZXXOLD(i)-ONE)/TWO
         ZYYOLD(i)=(ZYYOLD(i)-ONE)/TWO
         ZZZOLD(i)=(ZZZOLD(i)-ONE)/TWO
         ZXYOLD(i)=(ZXYOLD(i))/TWO
         ZYZOLD(i)=(ZYZOLD(i))/TWO
         ZZXOLD(i)=(ZZXOLD(i))/TWO
C
c        midpoint Green lagrange strains
C 
         ZXXMID(i)  =  (ZXX(i)+ ZXXOLD(i))/TWO
         ZYYMID(i)  =  (ZYY(i)+ ZYYOLD(i))/TWO
         ZZZMID(i)  =  (ZZZ(i)+ ZZZOLD(i))/TWO
         ZXYMID(i)  =  (ZXY(i)+ ZXYOLD(i))/TWO
         ZYZMID(i)  =  (ZYZ(i)+ ZYZOLD(i))/TWO
         ZZXMID(i)  =  (ZZX(i)+ ZZXOLD(i))/TWO
C
         JAC(I) = F(I,1,1)*F(I,2,2)*F(I,3,3) - F(I,1,1)*F(I,2,3)*F(I,3,2) -
     .          F(I,3,3)*F(I,1,2)*F(I,2,1)  + F(I,1,2)*F(I,2,3)*F(I,3,1) +
     .          F(I,2,1)*F(I,3,2)*F(I,1,3)  - F(I,2,2)*F(I,3,1)*F(I,1,3)
C
        print*, 'JAC(I) ', i, JAC(I)
       ENDDO
							
							
	!call conversion here in standalone

C
c     compute eigenvalues and eigenvectors
c     of the midpoint GL strain tensor
c     crucial : need to make sure that the order of eigenvalues
c     is the same as what is done in subroutine conversion
C

      DO I=1,NEL                       
        AV(1,I) = ZXXMID(I)
        AV(2,I) = ZYYMID(I)
        AV(3,I) = ZZZMID(I)
        AV(4,I) = ZXYMID(I)
        AV(5,I) = ZYZMID(I)
        AV(6,I) = ZZXMID(I)
      ENDDO                     
      IF (IRESP==1) THEN
          CALL VALPVECDP(AV,EVV,DIRPRV,NEL)
      ELSE
          CALL VALPVEC  (AV,EVV,DIRPRV,NEL)
      ENDIF

c     compute the strain rates in each one
c     of the principal directions of the midpoint GL tensor
c
c     just for the record : the eigenvectors of the midpoint GL
c     are of course not the eigenvectors of the strain rate tensor
c     this is a flaw of the formulation, kind of inevitable
c
c     the fields FSij are used here 
c     they will be used again in the stress transformation
c     these are just intermediate values that can be overwritten
c
c     VEC(1:3) is the first eigenvector
c     VEC(4:6) is the second eigenvector
c     VEC(7:9) is the third eigenvector
c     they may need to be reordered !!!
c
         !print*,'end Compute eigenvalues and eigenvector '
      DO I=1,NEL

          !calculation of D (strain rate must be checked- temporary)
          EP1 = EPSPXX(I)
          EP2 = EPSPYY(I)      
          EP3 = EPSPZZ(I) 
          EP4 = HALF*EPSPXY(I)        
          EP5 = HALF*EPSPYZ(I)
          EP6 = HALF*EPSPZX(I)
c c        phi_trans*L*phi_t    
c c        VEC-transposed * EPSP
          ERT11 =DIRPRV(1,1,I)*EP1 + DIRPRV(2,1,I)*EP4 + DIRPRV(3,1,I)*EP6
          ERT12 =DIRPRV(1,2,I)*EP1 + DIRPRV(2,2,I)*EP4 + DIRPRV(3,2,I)*EP6
          ERT13 =DIRPRV(1,3,I)*EP1 + DIRPRV(2,3,I)*EP4 + DIRPRV(3,3,I)*EP6
     
          ERT21 =DIRPRV(1,1,I)*EP4 + DIRPRV(2,1,I)*EP2 + DIRPRV(3,1,I)*EP5
          ERT22 =DIRPRV(1,2,I)*EP4 + DIRPRV(2,2,I)*EP2 + DIRPRV(3,2,I)*EP5
          ERT23 =DIRPRV(1,3,I)*EP4 + DIRPRV(2,3,I)*EP2 + DIRPRV(3,3,I)*EP5  
     
          ERT31 =DIRPRV(1,1,I)*EP6 + DIRPRV(2,1,I)*EP5 + DIRPRV(3,1,I)*EP3
          ERT32 =DIRPRV(1,2,I)*EP6 + DIRPRV(2,2,I)*EP5 + DIRPRV(3,2,I)*EP3
          ERT33 =DIRPRV(1,3,I)*EP6 + DIRPRV(2,3,I)*EP5 + DIRPRV(3,3,I)*EP3       
c C
c c       VEC-transposed * EPSP * VEC    diagonal components only
c C
          EPSP(1) = DIRPRV(1,1,I)*ERT11 + DIRPRV(2,1,I)*ERT21 
     .                                  + DIRPRV(3,1,I)*ERT31 
          EPSP(2) = DIRPRV(1,2,I)*ERT12 + DIRPRV(2,2,I)*ERT22 
     .                                  + DIRPRV(3,2,I)*ERT32 
          EPSP(3) = DIRPRV(1,3,I)*ERT13 + DIRPRV(2,3,I)*ERT23 
     .                                  + DIRPRV(3,3,I)*ERT33
           !(e > 0 compression and e < 0 traction)
                 print*, 'element ', i 

          print*, 'EVV(1,I) ,EVV(2,I) ', EVV(1,I) ,EVV(2,I)
          print*, '          EVV(3,I) ', EVV(3,I)
           EV(I,1) = SQRT(TWO * EVV(1,I) + ONE )
           EV(I,2) = SQRT(TWO * EVV(2,I) + ONE )
           EV(I,3) = SQRT(TWO * EVV(3,I) + ONE )

          STRAIN(I,1) = ONE - EV(i,1)        
          STRAIN(I,2) = ONE - EV(i,2)        
          STRAIN(I,3) = ONE - EV(i,3) 
c c       abs(eps) not necessary  
          DRATE(I,1) = EPSP(1)*(ONE - STRAIN(I,1)) ! eng
          DRATE(I,2) = EPSP(2)*(ONE - STRAIN(I,2))
          DRATE(I,3) = EPSP(3)*(ONE - STRAIN(I,3))
       ENDDO
c c
c c     compute global incremental stiffness matrix
c c     essentially a conversion from principal midpoint GL system
c c     to the global reference system
c c
c c     drate goes in to subroutine conversion as true strain rates
c c     in the principal directions of the midpoint GL strain tensor
c c
c     drate comes out of subroutine conversion as principal 2PK
c     dynamic overstress ( principal directions of GL and 2PK are same )
c
      CALL CONVERSION 
     &    (ZXXMID,ZYYMID,ZZZMID,ZXYMID,ZYZMID,ZZXMID,
     &       CIJKL,DFIRST,DRATE,
     &       NFUNC, IFUNC, SCAL, NEL, NPF, TF, JAC , SLOPE,
     &       NUMTABL ,TABLE)
C
      DO I = 1,NEL
C
C       increment of Green Lagrange strains
C
        DZXX(i)=ZXX(i)-ZXXOLD(i)
        DZYY(i)=ZYY(i)-ZYYOLD(i)
        DZZZ(i)=ZZZ(i)-ZZZOLD(i)
        DZXY(i)=ZXY(i)-ZXYOLD(i)
        DZYZ(i)=ZYZ(i)-ZYZOLD(i)
        DZZX(i)=ZZX(i)-ZZXOLD(i)
C
C       use Gamma's instead of epsilons for increments
C
        DZXY(i)=DZXY(i)*TWO
        DZYZ(i)=DZYZ(i)*TWO
        DZZX(i)=DZZX(i)*TWO
C
C       increment hyperelastic 2PK stresses   dS=C*dE
C
      DO 6 IK=1,6
        spknorate(IK,I)=spknorate(IK,I)
     &        + Cijkl(I,1,IK)*DZXX(I)+Cijkl(I,2,IK)*DZYY(I)
     &        + Cijkl(I,3,IK)*DZZZ(I)+Cijkl(I,4,IK)*DZXY(I)
     &        + Cijkl(I,5,IK)*DZYZ(I)+Cijkl(I,6,IK)*DZZX(I)
       !write(*,*) ' Cijkl(I,1,IK), Cijkl(I,2,IK) ', Cijkl(I,1,IK), Cijkl(I,2,IK) 
       !write(*,*) ' Cijkl(I,3,IK), Cijkl(I,4,IK) ', Cijkl(I,3,IK), Cijkl(I,4,IK)
       !write(*,*) ' Cijkl(I,5,IK), Cijkl(I,6,IK) ', Cijkl(I,5,IK), Cijkl(I,6,IK) 
 
       !write(*,*) '  spknorate(IK,I)',ik,  spknorate(IK,I)
   6  Continue 
       !write(*,*) 'DZXX(I), DZYY(I) ',DZXX(I), DZYY(I)
       !write(*,*) 'DZZZ(I), DZXY(I) ',DZZZ(I), DZXY(I) 
       !write(*,*) 'DZYZ(I), DZZX(I) ',DZYZ(I), DZZX(I)

          UVAR(I,11)  =  spknorate(1,i)         
          UVAR(I,12)  =  spknorate(2,i)         
          UVAR(I,13)  =  spknorate(3,i)         
          UVAR(I,14)  =  spknorate(4,i)         
          UVAR(I,15)  =  spknorate(5,i)         
          UVAR(I,16)  =  spknorate(6,i)         
C
c     add the viscous 2PK stresses ( dynamic overstress )
c     drate = dynamic principal 2PK - static principal 2PK
C

      do 7 IK=1,6
      if (IK.lt.4) then
      spknodam(IK,I)=spknorate(IK,I)
     &        +dfirst(1,IK,i)*drate(i,1)
     &        +dfirst(2,IK,i)*drate(i,2)
     &        +dfirst(3,IK,i)*drate(i,3)
      else
      spknodam(IK,I)=spknorate(IK,I)
     &        + (dfirst(1,IK,i)*drate(i,1)
     &        +  dfirst(2,IK,i)*drate(i,2)
     &        +  dfirst(3,IK,i)*drate(i,3) )*0.5d0
      endif
   7  continue 
c
      ENDDO
C
C     compute damage at t(n+1)
c    computes the global incremental stiffness matrix
c    at t(n+1/2) ( midpoint GL strains are used )

       CALL CONDAMAGE(ZXX ,ZYY  ,ZZZ ,ZXY      ,ZZX     ,ZYZ,
     .       DAMAGE ,NEL  ,HU  ,SHAPE   ,SCALINT,WHYSMAX,
     .       NUMTABL ,TABLE)
     
C
C     apply damage
C 
      DO I=1,NEL
        UVAR(I,17)    =  WHYSMAX(I)         
C    
        !print*, 'damage(i)', i, damage(i)
        do 8 IK=1,6
         !damage(i) = zero
          spk(IK,I)=spknodam(IK,I)*(ONE-damage(i))
   8    Continue
C
C       transform 2PK to Cauchy stresses  
C
        S1(i)=spk(1,i)                    
        S2(i)=spk(2,i)                    
        S3(i)=spk(3,i)                    
        S4(i)=spk(4,i)                    
        S5(i)=spk(5,i)                    
        S6(i)=spk(6,i)                    
c
        FS11(i) = F(I,1,1)*S1(i) + F(I,1,2)*S4(i) + F(I,1,3)*S6(i)
        FS12(i) = F(I,1,1)*S4(i) + F(I,1,2)*S2(i) + F(I,1,3)*S5(i)
        FS13(i) = F(I,1,1)*S6(i) + F(I,1,2)*S5(i) + F(I,1,3)*S3(i)
        FS21(i) = F(I,2,1)*S1(i) + F(I,2,2)*S4(i) + F(I,2,3)*S6(i)
        FS22(i) = F(I,2,1)*S4(i) + F(I,2,2)*S2(i) + F(I,2,3)*S5(i)
        FS23(i) = F(I,2,1)*S6(i) + F(I,2,2)*S5(i) + F(I,2,3)*S3(i)
        FS31(i) = F(I,3,1)*S1(i) + F(I,3,2)*S4(i) + F(I,3,3)*S6(i)
        FS32(i) = F(I,3,1)*S4(i) + F(I,3,2)*S2(i) + F(I,3,3)*S5(i)
        FS33(i) = F(I,3,1)*S6(i) + F(I,3,2)*S5(i) + F(I,3,3)*S3(i)
c
        SIGNXX(I) = FS11(i)*F(I,1,1)+FS12(i)*F(I,1,2)+FS13(i)*F(I,1,3)
        SIGNYY(I) = FS21(i)*F(I,2,1)+FS22(i)*F(I,2,2)+FS23(i)*F(I,2,3)
        SIGNZZ(I) = FS31(i)*F(I,3,1)+FS32(i)*F(I,3,2)+FS33(i)*F(I,3,3)
        SIGNXY(I) = FS11(i)*F(I,2,1)+FS12(i)*F(I,2,2)+FS13(i)*F(I,2,3)
        SIGNYZ(I) = FS21(i)*F(I,3,1)+FS22(i)*F(I,3,2)+FS23(i)*F(I,3,3)
        SIGNZX(I) = FS11(i)*F(I,3,1)+FS12(i)*F(I,3,2)+FS13(i)*F(I,3,3)
C
        SIGNXX(I) = SIGNXX(I) / JAC(I)   
        SIGNYY(I) = SIGNYY(I) / JAC(I)   
        SIGNZZ(I) = SIGNZZ(I) / JAC(I)   
        SIGNXY(I) = SIGNXY(I) / JAC(I)   
        SIGNYZ(I) = SIGNYZ(I) / JAC(I)   
        SIGNZX(I) = SIGNZX(I) / JAC(I)  
        !print*, 'damage(i)', i, damage(i)
        !write(*,*) 'SIGNXX, SIGNYY ,SIGNZZ ' , SIGNXX(I), SIGNYY(I) ,SIGNZZ(I)

C
C       timestep if the slope of the stress/strain curve
C       exceeds the user specified modulus
C
         GS=G
         KS=K
         SLOPEMAX(I)=MAX(SLOPE(I,1), SLOPE(I,2))
         SLOPEMAX(I)=MAX(SLOPEMAX(I), SLOPE(I,3))
         IF (SLOPEMAX(I).GT.EMOD) THEN
           GS=SLOPEMAX(I)/2.
           KS=SLOPEMAX(I)/3.
         ENDIF
         SOUNDSP(I) = SQRT((4.0*GS/3.0 + KS)/RHO(I))
         VISCMAX(I) = 0.
c        print *, I ,' soundsp '
      ENDDO

C Output data structure
C-------------------------------------------
C 
C-------------------------------------------
C
C-------------------------------------------------------------------------
      RETURN
      END
