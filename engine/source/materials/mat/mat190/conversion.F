c===========================================================================
       SUBROUTINE CONVERSION 
     .    (ZXX,ZYY,ZZZ,ZXY,ZYZ,ZZX,CIJKL,
     .      DFIRST,DRATE,
     .      NFUNC,IFUNC, SCAL, NEL, NPF, TF, JAC, SLOPE ,NUMTABL ,TABLE)
c           computation of coefficients to convert stiffness
c           matrix from local to global reference system
c===========================================================================
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE TABLE_MOD
      USE INTERFACE_TABLE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "comlock.inc"
C-----------------------------------------------
#include      "scr05_c.inc"
#include      "impl1_c.inc"
#include      "tabsiz_c.inc"
c===========================================================================
C-----------------------------------------------
C   I N T P U T   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) ::NEL,   NFUNC ,NUMTABL
      my_real, INTENT(IN) :: ZXX(NEL),ZYY(NEL),ZZZ(NEL),ZXY(NEL),ZZX(NEL),ZYZ(NEL), JAC(NEL)
      INTEGER ,DIMENSION(NFUNC),INTENT(IN) :: IFUNC
      INTEGER ,DIMENSION(SNPC) ,INTENT(IN) :: NPF
      my_real ,DIMENSION(STF)  ,INTENT(IN) :: TF
C-----------------------------------------------
C    I N T P U T   O U T P U T   A r g u m e n t s
C-----------------------------------------------
      my_real, INTENT(OUT) :: DFIRST(3,6,NEL),CIJKL(NEL,6,6), SLOPE(NEL,3)
      my_real, INTENT(INOUT) :: DRATE(NEL,3)
      my_real, EXTERNAL :: FINTER
      TYPE (TTABLE), DIMENSION(NUMTABL) ,TARGET  ::  TABLE
c===========================================================================
C   EXTERNAL FINTER
C        Y = FINTER(IFUNC(J),X,NPF,TF,DYDX)
C        Y       : y = f(x)
C        X       : x
C        DYDX    : f'(x) = dy/dx
C        IFUNC(J): FUNCTION INDEX
C              J : FIRST(J=1), SECOND(J=2) .. FUNCTION USED FOR THIS LAW
C        NPF,TF  : FUNCTION PARAMETER
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  I,J,K,L
 
      my_real, DIMENSION(NEL)   :: DYDX1,VOLSTR, !TRUERATE, 
     .                             FSIG1,DSIG1,FSIG2,DSIG2,
     .                             SIGENER,DSIGENER,FSIG3,DSIG3
      my_real, DIMENSION(NEL,3)   :: ALAMBDA,EPSILON,ENGIRATE, SIGMA, DYNSIGMA,
     .                              STRAINRATE 
      my_real
     .   SCAL, 
     .    DIK(NEL,3), SIK(NEL,3), 
     .    ac,as,Z1(NEL),Z2(NEL),Z3(NEL),dsecond (3,6,6,NEL),
     .   ah1(NEL),ah2(NEL),ah3(NEL),at2(NEL),at3(NEL),aa3(NEL),aa2(NEL),aa4(NEL),
     .   dzero(3,NEL),
     .   dat2(6,NEL),
     .   dat3(6,NEL),
     .   da(6,NEL),
     .   da2(6,6,NEL),
     .   dat22(6,6,NEL),
     .   dat32(6,6,NEL),
     .   dx2(6,6,NEL),
     .  dh1Zxx,dh1Zyy,dh1Zxy,dh1Zyz,dh1Zzx,dh1Zzz,
     .  dh2Zxx,dh2Zyy,dh2Zxy,dh2Zyz,dh2Zzx,dh2Zzz,
     .  dh3Zxx,dh3Zyy,dh3Zxy,dh3Zyz,dh3Zzx,dh3Zzz,
     .  dt2Zxx,dt2Zyy,dt2Zxy,dt2Zyz,dt2Zzx,dt2Zzz,
     .  dt3Zxx,dt3Zyy,dt3Zxy,dt3Zyz,dt3Zzx,dt3Zzz,
     .  temp1,temp2,temp3,temp4,temp5,temp11,temp12,
     .  dZ1Zxx,dZ1Zyy,dZ1Zxy,dZ1Zyz,dZ1Zzx,dZ1Zzz,
     .  dZ2Zxx,dZ2Zyy,dZ2Zxy,dZ2Zyz,dZ2Zzx,dZ2Zzz,
     .  dZ3Zxx,dZ3Zyy,dZ3Zxy,dZ3Zyz,dZ3Zzx,dZ3Zzz,
     .  da3Zxx,da3Zyy,da3Zxy,da3Zyz,da3Zzx,da3Zzz,
     .   x1,x2,x22,x3,x4,x5,x6,
     .   t1,t2,t3,t4,t5,t6,
     .   afac1,afac2,afac3,afac4,afac33,afac5,
     .   ahelp3,ahelp33,ahelp,ahelp4,
     .   xcap,dxcap,acosteta,
     .   dah2ZxxZxx,dah2ZxxZyy,dah2ZxxZzz,dah2ZxxZxy,dah2ZxxZyz,dah2ZxxZzx,
     .   dah2ZyyZxx,dah2ZyyZyy,dah2ZyyZzz,dah2ZyyZxy,dah2ZyyZyz,dah2ZyyZzx,
     .   dah2ZzzZxx,dah2ZzzZyy,dah2ZzzZzz,dah2ZzzZxy,dah2ZzzZyz,dah2ZzzZzx,
     .   dah2ZxyZxx,dah2ZxyZyy,dah2ZxyZzz,dah2ZxyZxy,dah2ZxyZyz,dah2ZxyZzx,
     .   dah2ZyzZxx,dah2ZyzZyy,dah2ZyzZzz,dah2ZyzZxy,dah2ZyzZyz,dah2ZyzZzx,
     .   dah2ZzxZxx,dah2ZzxZyy,dah2ZzxZzz,dah2ZzxZxy,dah2ZzxZyz,dah2ZzxZzx,
     .   dah3ZxxZxx,dah3ZxxZyy,dah3ZxxZzz,dah3ZxxZxy,dah3ZxxZyz,dah3ZxxZzx,
     .   dah3ZyyZxx,dah3ZyyZyy,dah3ZyyZzz,dah3ZyyZxy,dah3ZyyZyz,dah3ZyyZzx,
     .   dah3ZzzZxx,dah3ZzzZyy,dah3ZzzZzz,dah3ZzzZxy,dah3ZzzZyz,dah3ZzzZzx,
     .   dah3ZxyZxx,dah3ZxyZyy,dah3ZxyZzz,dah3ZxyZxy,dah3ZxyZyz,dah3ZxyZzx,
     .   dah3ZyzZxx,dah3ZyzZyy,dah3ZyzZzz,dah3ZyzZxy,dah3ZyzZyz,dah3ZyzZzx,
     .   dah3ZzxZxx,dah3ZzxZyy,dah3ZzxZzz,dah3ZzxZxy,dah3ZzxZyz,dah3ZzxZzx
      ! TABLE
      INTEGER  IPOS(NEL,3)
      my_real, DIMENSION(NEL,3)   :: XVEC
      TYPE(TTABLE), POINTER :: FUNC_ENER , FUNC_SIG
c===========================================================================
c===========================================================================
      FUNC_SIG   => TABLE(1)
      FUNC_ENER  => TABLE(2)
      

      do i=1,NEL
c      
c     part 1 : principal strains
c
c     compute H1, H2 and H3
c
      ah1(i)=Zxx(i)+Zyy(i)+Zzz(i)
      ah2(i)=Zyy(i)*Zzz(i)+Zzz(i)*Zxx(i)+Zxx(i)*Zyy(i)
     .        -Zyz(i)**2-Zzx(i)**2-Zxy(i)**2
      ah3(i)=Zxx(i)*Zyy(i)*Zzz(i)+2.d0*Zxy(i)*Zyz(i)*Zzx(i)
     .        -Zxx(i)*Zyz(i)**2-Zyy(i)*Zzx(i)**2-Zzz(i)*Zxy(i)**2
c
c    derivatives of h1
c
      dh1Zxx=1.d0
      dh1Zyy=1.d0
      dh1Zzz=1.d0
      dh1Zxy=0.
      dh1Zyz=0.
      dh1Zzx=0.
c
c    derivatives of h2
c
      dh2Zxx=Zyy(i)+Zzz(i)
      dh2Zyy=Zxx(i)+Zzz(i)
      dh2Zzz=Zyy(i)+Zxx(i)
      dh2Zxy=-2.d0*Zxy(i)
      dh2Zyz=-2.d0*Zyz(i)
      dh2Zzx=-2.d0*Zzx(i)
c
c    derivatives of h3
c
      dh3Zxx=Zyy(i)*Zzz(i)-Zyz(i)**2
      dh3Zyy=Zxx(i)*Zzz(i)-Zzx(i)**2
      dh3Zzz=Zyy(i)*Zxx(i)-Zxy(i)**2
      dh3Zxy=-2.d0*Zxy(i)*Zzz(i)+2.d0*Zyz(i)*Zzx(i)
      dh3Zyz=-2.d0*Zyz(i)*Zxx(i)+2.d0*Zxy(i)*Zzx(i)
      dh3Zzx=-2.d0*Zzx(i)*Zyy(i)+2.d0*Zyz(i)*Zxy(i)
c
c     compute T2 and T3
c
      at3(i)=2.d0*ah1(i)**3/TWENTY7 - ah1(i)*ah2(i)/3.d0 +  ah3(i)
      at2(i)=-ah2(i) + ah1(i)**2/3.d0
      IF(at2(i) <= zero) at2(i)=em20 !max(at2(i),1.e-16)
c
c     compute dcos(a)
c
      aa3(i)=at3(i)/2.d0/dsqrt(at2(i)**3/27.d0)
      !print*, 'aa3(i) ',aa3(i)  
      aa3(i)=min(1.0d0,aa3(i)) 
      aa3(i)=max(-1.0d0,aa3(i)) 
      aa2(i)=dsqrt(at2(i)/3.d0)
      acosteta = dacos(aa3(i))
c
c     compute principal strains
c
      pi=3.14159265359d0
      Z1(i)=TWO*aa2(i)*dcos((THIRD) * acosteta  ) 
      Z2(i)=TWO*aa2(i)*dcos((THIRD) *(acosteta -TWO *pi) ) 
      Z3(i)=TWO*aa2(i)*dcos((THIRD) *(acosteta -FOUR*pi) ) 
      !print*, 'element ', i 
      !print*, 'Z1(i) , Z2(i) ',Z1(i) , Z2(i) 
      !print*, 'Z3(i)         ',Z3(i)
c
c     add hydrostatic part
c
      Z1(i)=Z1(i)+ah1(i)/3.d0
      Z2(i)=Z2(i)+ah1(i)/3.d0
      Z3(i)=Z3(i)+ah1(i)/3.d0  
      !print*, 'hydro ', i 
      !print*, 'Z1(i) , Z2(i) ',Z1(i) , Z2(i) 
      !print*, 'Z3(i)         ',Z3(i)
c
c
c
c     part 2
c     first derivatives of principal strains
c
c
c     derivatives  of  t2
c
      dt2Zxx=2.d0*ah1(i)/3.d0-(Zyy(i)+Zzz(i))
      dt2Zyy=2.d0*ah1(i)/3.d0-(Zxx(i)+Zzz(i))
      dt2Zzz=2.d0*ah1(i)/3.d0-(Zyy(i)+Zxx(i))
      dt2Zxy=2.d0*Zxy(i)
      dt2Zyz=2.d0*Zyz(i)
      dt2Zzx=2.d0*Zzx(i)
c
c     store derivatives of t2
c
      dat2(1,i)=dt2Zxx
      dat2(2,i)=dt2Zyy
      dat2(3,i)=dt2Zzz
      dat2(4,i)=dt2Zxy
      dat2(5,i)=dt2Zyz
      dat2(6,i)=dt2Zzx
c
c     derivatives of t3
c    
      dt3Zxx=6.d0*ah1(i)**2/27.d0-ah2(i)/3.d0
     .                 -ah1(i)*(Zyy(i)+Zzz(i))/3.d0
     .                 +(Zyy(i)*Zzz(i)-Zyz(i)**2)
      dt3Zyy=6.d0*ah1(i)**2/27.d0-ah2(i)/3.d0
     .                 -ah1(i)*(Zxx(i)+Zzz(i))/3.d0
     .                 +(Zxx(i)*Zzz(i)-Zzx(i)**2)
      dt3Zzz=6.d0*ah1(i)**2/27.d0-ah2(i)/3.d0
     .                 -ah1(i)*(Zyy(i)+Zxx(i))/3.d0
     .                 +(Zyy(i)*Zxx(i)-Zxy(i)**2)
      dt3Zxy=2.d0*ah1(i)*Zxy(i)/3.d0
     .       +2.d0*(Zyz(i)*Zzx(i)-Zxy(i)*Zzz(i))
      dt3Zyz=2.d0*ah1(i)*Zyz(i)/3.d0
     .       +2.d0*(Zxy(i)*Zzx(i)-Zyz(i)*Zxx(i))
      dt3Zzx=2.d0*ah1(i)*Zzx(i)/3.d0
     .       +2.d0*(Zyz(i)*Zxy(i)-Zzx(i)*Zyy(i))
c
c     store derivatives of t3
c
      dat3(1,i)=dt3Zxx
      dat3(2,i)=dt3Zyy
      dat3(3,i)=dt3Zzz
      dat3(4,i)=dt3Zxy
      dat3(5,i)=dt3Zyz
      dat3(6,i)=dt3Zzx
c
c     auxillary coefficients
c
      temp11=1.d0-(at3(i)/2.d0/dsqrt(at2(i)**3/27.d0))**2
      temp12=max(temp11,1.e-16)
      temp1=-1.d0/dsqrt(temp12)
c
c     'if' test needed to get correct limit value
c     in the uniaxial cases
c     biaxial cases seem fine without the test
c
      if (temp11.lt.1.e-16) temp1=0.0d0
c
      temp2=1.d0/2.d0/dsqrt(at2(i)**3/27.d0)
      temp3=(-1.d0/2.d0)*at3(i)*3.d0*at2(i)**2/2.d0/27.d0
     .      /dsqrt((at2(i)**3/27.d0)**3)
c
c     derivatives of a
c
      da3Zxx=temp1*(temp2*dt3Zxx+temp3*dt2Zxx) 
      da3Zyy=temp1*(temp2*dt3Zyy+temp3*dt2Zyy) 
      da3Zzz=temp1*(temp2*dt3Zzz+temp3*dt2Zzz) 
      da3Zxy=temp1*(temp2*dt3Zxy+temp3*dt2Zxy) 
      da3Zyz=temp1*(temp2*dt3Zyz+temp3*dt2Zyz) 
      da3Zzx=temp1*(temp2*dt3Zzx+temp3*dt2Zzx) 
c
c     store derivatives of a
c
      da(1,i)=da3Zxx
      da(2,i)=da3Zyy
      da(3,i)=da3Zzz
      da(4,i)=da3Zxy
      da(5,i)=da3Zyz
      da(6,i)=da3Zzx
c
      temp4=1.d0/3.d0/dsqrt(at2(i)/3.d0)
      temp5=2.d0*dsqrt(at2(i)/3.d0)/3.d0
      aa4(i)=dacos(aa3(i))
c
c     derivatives of principal strains
c
      dZ1Zxx=temp4*dt2Zxx*dcos(-aa4(i)/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0)*da3Zxx
      dZ2Zxx=temp4*dt2Zxx*dcos(-aa4(i)/3.d0+2*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+2*pi/3.d0)*da3Zxx
      dZ3Zxx=temp4*dt2Zxx*dcos(-aa4(i)/3.+4*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+4*pi/3.d0)*da3Zxx
      dZ1Zyy=temp4*dt2Zyy*dcos(-aa4(i)/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0)*da3Zyy
      dZ2Zyy=temp4*dt2Zyy*dcos(-aa4(i)/3.d0+2*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+2*pi/3.d0)*da3Zyy
      dZ3Zyy=temp4*dt2Zyy*dcos(-aa4(i)/3.d0+4*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+4*pi/3.d0)*da3Zyy
      dZ1Zzz=temp4*dt2Zzz*dcos(-aa4(i)/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0)*da3Zzz
      dZ2Zzz=temp4*dt2Zzz*dcos(-aa4(i)/3.d0+2*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+2*pi/3.d0)*da3Zzz
      dZ3Zzz=temp4*dt2Zzz*dcos(-aa4(i)/3.d0+4*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+4*pi/3.d0)*da3Zzz
      dZ1Zxy=temp4*dt2Zxy*dcos(-aa4(i)/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0)*da3Zxy
      dZ2Zxy=temp4*dt2Zxy*dcos(-aa4(i)/3.d0+2*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+2*pi/3.d0)*da3Zxy
      dZ3Zxy=temp4*dt2Zxy*dcos(-aa4(i)/3.d0+4*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+4*pi/3.d0)*da3Zxy
      dZ1Zyz=temp4*dt2Zyz*dcos(-aa4(i)/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0)*da3Zyz
      dZ2Zyz=temp4*dt2Zyz*dcos(-aa4(i)/3.d0+2*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+2*pi/3.d0)*da3Zyz
      dZ3Zyz=temp4*dt2Zyz*dcos(-aa4(i)/3.d0+4*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+4*pi/3.d0)*da3Zyz
      dZ1Zzx=temp4*dt2Zzx*dcos(-aa4(i)/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0)*da3Zzx
      dZ2Zzx=temp4*dt2Zzx*dcos(-aa4(i)/3.d0+2*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+2*pi/3.d0)*da3Zzx
      dZ3Zzx=temp4*dt2Zzx*dcos(-aa4(i)/3.d0+4*pi/3.d0)
     .                 +temp5*dsin(-aa4(i)/3.d0+4*pi/3.d0)*da3Zzx
c
c     store first derivatives
c     store zero'd derivatives
c
      dzero(1,i)=Z1(i)
      dzero(2,i)=Z2(i)
      dzero(3,i)=Z3(i)
      dfirst(1,1,i)=dZ1Zxx
      dfirst(1,2,i)=dZ1Zyy
      dfirst(1,3,i)=dZ1Zzz
      dfirst(1,4,i)=dZ1Zxy
      dfirst(1,5,i)=dZ1Zyz
      dfirst(1,6,i)=dZ1Zzx
      dfirst(2,1,i)=dZ2Zxx
      dfirst(2,2,i)=dZ2Zyy
      dfirst(2,3,i)=dZ2Zzz
      dfirst(2,4,i)=dZ2Zxy
      dfirst(2,5,i)=dZ2Zyz
      dfirst(2,6,i)=dZ2Zzx
      dfirst(3,1,i)=dZ3Zxx
      dfirst(3,2,i)=dZ3Zyy
      dfirst(3,3,i)=dZ3Zzz
      dfirst(3,4,i)=dZ3Zxy
      dfirst(3,5,i)=dZ3Zyz
      dfirst(3,6,i)=dZ3Zzx
c
c    add hydrostatic part
c
      dfirst(1,1,i)=dZ1Zxx+1.d0/3.
      dfirst(1,2,i)=dZ1Zyy+1.d0/3.
      dfirst(1,3,i)=dZ1Zzz+1.d0/3.
      dfirst(2,1,i)=dZ2Zxx+1.d0/3.
      dfirst(2,2,i)=dZ2Zyy+1.d0/3.
      dfirst(2,3,i)=dZ2Zzz+1.d0/3.
      dfirst(3,1,i)=dZ3Zxx+1.d0/3.
      dfirst(3,2,i)=dZ3Zyy+1.d0/3.
      dfirst(3,3,i)=dZ3Zzz+1.d0/3.
c
c     code as in MAT_133 up till here
c     so this part I assume works
c
c
c     part 3
c     second derivatives of principal strains
c
c
c     second derivatives of H1 are all zero
c
c
c     second derivatives of H2
c
      dah2ZxxZxx=0.
      dah2ZxxZyy=1.d0
      dah2ZxxZzz=1.d0
      dah2ZxxZxy=0.
      dah2ZxxZyz=0.
      dah2ZxxZzx=0.
c
      dah2ZyyZxx=1.d0
      dah2ZyyZyy=0.
      dah2ZyyZzz=1.d0
      dah2ZyyZxy=0.
      dah2ZyyZyz=0.
      dah2ZyyZzx=0.
c
      dah2ZzzZxx=1.d0
      dah2ZzzZyy=1.d0
      dah2ZzzZzz=0.
      dah2ZzzZxy=0.
      dah2ZzzZyz=0.
      dah2ZzzZzx=0.
c
      dah2ZxyZxx=0.
      dah2ZxyZyy=0.
      dah2ZxyZzz=0.
      dah2ZxyZxy=-2.d0
      dah2ZxyZyz=0.
      dah2ZxyZzx=0.
c
      dah2ZyzZxx=0.
      dah2ZyzZyy=0.
      dah2ZyzZzz=0.
      dah2ZyzZxy=0.
      dah2ZyzZyz=-2.d0
      dah2ZyzZzx=0.
c
      dah2ZzxZxx=0.
      dah2ZzxZyy=0.
      dah2ZzxZzz=0.
      dah2ZzxZxy=0.
      dah2ZzxZyz=0.
      dah2ZzxZzx=-2.d0
c
c     second derivatives of H3
c
      dah3ZxxZxx=0.d0
      dah3ZxxZyy=Zzz(i)
      dah3ZxxZzz=Zyy(i)
      dah3ZxxZxy=0.d0
      dah3ZxxZyz=-2.d0*Zyz(i)
      dah3ZxxZzx=0.d0
c
      dah3ZyyZxx=Zzz(i)
      dah3ZyyZyy=0.d0
      dah3ZyyZzz=Zxx(i)
      dah3ZyyZxy=0.d0
      dah3ZyyZyz=0.d0
      dah3ZyyZzx=-2.d0*Zzx(i)
c
      dah3ZzzZxx=Zyy(i)
      dah3ZzzZyy=Zxx(i)
      dah3ZzzZzz=0.d0
      dah3ZzzZxy=-2.d0*Zxy(i)
      dah3ZzzZyz=0.d0
      dah3ZzzZzx=0.d0
c
      dah3ZxyZxx=0.d0
      dah3ZxyZyy=0.d0
      dah3ZxyZzz=-2.d0*Zxy(i)
      dah3ZxyZxy=-2.d0*Zzz(i)
      dah3ZxyZyz=2.d0*Zzx(i)
      dah3ZxyZzx=2.d0*Zyz(i)
c
      dah3ZyzZxx=-2.d0*Zyz(i)
      dah3ZyzZyy=0.d0
      dah3ZyzZzz=0.d0
      dah3ZyzZxy=2.d0*Zzx(i)
      dah3ZyzZyz=-2.d0*Zxx(i)
      dah3ZyzZzx=2.d0*Zxy(i)
c
      dah3ZzxZxx=0.d0
      dah3ZzxZyy=-2.d0*Zzx(i)
      dah3ZzxZzz=0.d0
      dah3ZzxZxy=2.d0*Zyz(i)
      dah3ZzxZyz=2.d0*Zxy(i)
      dah3ZzxZzx=-2.d0*Zyy(i)
c
c     second derivatives of T2
c
      dat22(1,1,i)=-dah2ZxxZxx+(2.d0/3.d0)
      dat22(1,2,i)=-dah2ZxxZyy+(2.d0/3.d0)
      dat22(1,3,i)=-dah2ZxxZzz+(2.d0/3.d0)
      dat22(1,4,i)=-dah2ZxxZxy
      dat22(1,5,i)=-dah2ZxxZyz
      dat22(1,6,i)=-dah2ZxxZzx
c
      dat22(2,1,i)=-dah2ZyyZxx+(2.d0/3.d0)
      dat22(2,2,i)=-dah2ZyyZyy+(2.d0/3.d0)
      dat22(2,3,i)=-dah2ZyyZzz+(2.d0/3.d0)
      dat22(2,4,i)=-dah2ZyyZxy
      dat22(2,5,i)=-dah2ZyyZyz
      dat22(2,6,i)=-dah2ZyyZzx
c
      dat22(3,1,i)=-dah2ZzzZxx+(2.d0/3.d0)
      dat22(3,2,i)=-dah2ZzzZyy+(2.d0/3.d0)
      dat22(3,3,i)=-dah2ZzzZzz+(2.d0/3.d0)
      dat22(3,4,i)=-dah2ZzzZxy
      dat22(3,5,i)=-dah2ZzzZyz
      dat22(3,6,i)=-dah2ZzzZzx
c
      dat22(4,1,i)=-dah2ZxyZxx
      dat22(4,2,i)=-dah2ZxyZyy
      dat22(4,3,i)=-dah2ZxyZzz
      dat22(4,4,i)=-dah2ZxyZxy
      dat22(4,5,i)=-dah2ZxyZyz
      dat22(4,6,i)=-dah2ZxyZzx
c
      dat22(5,1,i)=-dah2ZyzZxx
      dat22(5,2,i)=-dah2ZyzZyy
      dat22(5,3,i)=-dah2ZyzZzz
      dat22(5,4,i)=-dah2ZyzZxy
      dat22(5,5,i)=-dah2ZyzZyz
      dat22(5,6,i)=-dah2ZyzZzx
c
      dat22(6,1,i)=-dah2ZzxZxx
      dat22(6,2,i)=-dah2ZzxZyy
      dat22(6,3,i)=-dah2ZzxZzz
      dat22(6,4,i)=-dah2ZzxZxy
      dat22(6,5,i)=-dah2ZzxZyz
      dat22(6,6,i)=-dah2ZzxZzx
c
c     second derivatives of T3
c
      dat32(1,1,i)=dah3ZxxZxx-(ah1(i)/3.d0)*dah2ZxxZxx
      dat32(1,2,i)=dah3ZxxZyy-(ah1(i)/3.d0)*dah2ZxxZyy
      dat32(1,3,i)=dah3ZxxZzz-(ah1(i)/3.d0)*dah2ZxxZzz
      dat32(1,4,i)=dah3ZxxZxy-(ah1(i)/3.d0)*dah2ZxxZxy
      dat32(1,5,i)=dah3ZxxZyz-(ah1(i)/3.d0)*dah2ZxxZyz
      dat32(1,6,i)=dah3ZxxZzx-(ah1(i)/3.d0)*dah2ZxxZzx
c
      dat32(2,1,i)=dah3ZyyZxx-(ah1(i)/3.d0)*dah2ZyyZxx
      dat32(2,2,i)=dah3ZyyZyy-(ah1(i)/3.d0)*dah2ZyyZyy
      dat32(2,3,i)=dah3ZyyZzz-(ah1(i)/3.d0)*dah2ZyyZzz
      dat32(2,4,i)=dah3ZyyZxy-(ah1(i)/3.d0)*dah2ZyyZxy
      dat32(2,5,i)=dah3ZyyZyz-(ah1(i)/3.d0)*dah2ZyyZyz
      dat32(2,6,i)=dah3ZyyZzx-(ah1(i)/3.d0)*dah2ZyyZzx
c
      dat32(3,1,i)=dah3ZzzZxx-(ah1(i)/3.d0)*dah2ZzzZxx
      dat32(3,2,i)=dah3ZzzZyy-(ah1(i)/3.d0)*dah2ZzzZyy
      dat32(3,3,i)=dah3ZzzZzz-(ah1(i)/3.d0)*dah2ZzzZzz
      dat32(3,4,i)=dah3ZzzZxy-(ah1(i)/3.d0)*dah2ZzzZxy
      dat32(3,5,i)=dah3ZzzZyz-(ah1(i)/3.d0)*dah2ZzzZyz
      dat32(3,6,i)=dah3ZzzZzx-(ah1(i)/3.d0)*dah2ZzzZzx
c
      dat32(4,1,i)=dah3ZxyZxx-(ah1(i)/3.d0)*dah2ZxyZxx
      dat32(4,2,i)=dah3ZxyZyy-(ah1(i)/3.d0)*dah2ZxyZyy
      dat32(4,3,i)=dah3ZxyZzz-(ah1(i)/3.d0)*dah2ZxyZzz
      dat32(4,4,i)=dah3ZxyZxy-(ah1(i)/3.d0)*dah2ZxyZxy
      dat32(4,5,i)=dah3ZxyZyz-(ah1(i)/3.d0)*dah2ZxyZyz
      dat32(4,6,i)=dah3ZxyZzx-(ah1(i)/3.d0)*dah2ZxyZzx
c
      dat32(5,1,i)=dah3ZyzZxx-(ah1(i)/3.d0)*dah2ZyzZxx
      dat32(5,2,i)=dah3ZyzZyy-(ah1(i)/3.d0)*dah2ZyzZyy
      dat32(5,3,i)=dah3ZyzZzz-(ah1(i)/3.d0)*dah2ZyzZzz
      dat32(5,4,i)=dah3ZyzZxy-(ah1(i)/3.d0)*dah2ZyzZxy
      dat32(5,5,i)=dah3ZyzZyz-(ah1(i)/3.d0)*dah2ZyzZyz
      dat32(5,6,i)=dah3ZyzZzx-(ah1(i)/3.d0)*dah2ZyzZzx

      dat32(6,1,i)=dah3ZzxZxx-(ah1(i)/3.d0)*dah2ZzxZxx
      dat32(6,2,i)=dah3ZzxZyy-(ah1(i)/3.d0)*dah2ZzxZyy
      dat32(6,3,i)=dah3ZzxZzz-(ah1(i)/3.d0)*dah2ZzxZzz
      dat32(6,4,i)=dah3ZzxZxy-(ah1(i)/3.d0)*dah2ZzxZxy
      dat32(6,5,i)=dah3ZzxZyz-(ah1(i)/3.d0)*dah2ZzxZyz
      dat32(6,6,i)=dah3ZzxZzx-(ah1(i)/3.d0)*dah2ZzxZzx
c
c     add first derivative terms
c
      dat32(1,1,i)=dat32(1,1,i)+(12.d0*ah1(i)/27.d0)*dh1Zxx*dh1Zxx
     .       -(1.d0/3.d0)*dh1Zxx*dh2Zxx-(1.d0/3.d0)*dh1Zxx*dh2Zxx
      dat32(1,2,i)=dat32(1,2,i)+(12.d0*ah1(i)/27.d0)*dh1Zyy*dh1Zxx
     .       -(1.d0/3.d0)*dh1Zxx*dh2Zyy-(1.d0/3.d0)*dh1Zyy*dh2Zxx
      dat32(1,3,i)=dat32(1,3,i)+(12.d0*ah1(i)/27.d0)*dh1Zzz*dh1Zxx
     .       -(1.d0/3.d0)*dh1Zxx*dh2Zzz-(1.d0/3.d0)*dh1Zzz*dh2Zxx
      dat32(1,4,i)=dat32(1,4,i)+(12.d0*ah1(i)/27.d0)*dh1Zxy*dh1Zxx
     .       -(1.d0/3.d0)*dh1Zxx*dh2Zxy-(1.d0/3.d0)*dh1Zxy*dh2Zxx
      dat32(1,5,i)=dat32(1,5,i)+(12.d0*ah1(i)/27.d0)*dh1Zyz*dh1Zxx
     .       -(1.d0/3.d0)*dh1Zxx*dh2Zyz-(1.d0/3.d0)*dh1Zyz*dh2Zxx
      dat32(1,6,i)=dat32(1,6,i)+(12.d0*ah1(i)/27.d0)*dh1Zzx*dh1zxx
     .       -(1.d0/3.d0)*dh1Zxx*dh2Zzx-(1.d0/3.d0)*dh1Zzx*dh2Zxx
c
      dat32(2,1,i)=dat32(2,1,i)+(12.d0*ah1(i)/27.d0)*dh1Zxx*dh1Zyy
     .       -(1.d0/3.d0)*dh1Zyy*dh2Zxx-(1.d0/3.d0)*dh1Zxx*dh2Zyy
      dat32(2,2,i)=dat32(2,2,i)+(12.d0*ah1(i)/27.d0)*dh1Zyy*dh1Zyy
     .       -(1.d0/3.d0)*dh1Zyy*dh2Zyy-(1.d0/3.d0)*dh1Zyy*dh2Zyy
      dat32(2,3,i)=dat32(2,3,i)+(12.d0*ah1(i)/27.d0)*dh1Zzz*dh1Zyy
     .       -(1.d0/3.d0)*dh1Zyy*dh2Zzz-(1.d0/3.d0)*dh1Zzz*dh2Zyy
      dat32(2,4,i)=dat32(2,4,i)+(12.d0*ah1(i)/27.d0)*dh1Zxy*dh1Zyy
     .       -(1.d0/3.d0)*dh1Zyy*dh2Zxy-(1.d0/3.d0)*dh1Zxy*dh2Zyy
      dat32(2,5,i)=dat32(2,5,i)+(12.d0*ah1(i)/27.d0)*dh1Zyz*dh1Zyy
     .       -(1.d0/3.d0)*dh1Zyy*dh2Zyz-(1.d0/3.d0)*dh1Zyz*dh2Zyy
      dat32(2,6,i)=dat32(2,6,i)+(12.d0*ah1(i)/27.d0)*dh1Zzx*dh1zyy
     .       -(1.d0/3.d0)*dh1Zyy*dh2Zzx-(1.d0/3.d0)*dh1Zzx*dh2Zyy
c
      dat32(3,1,i)=dat32(3,1,i)+(12.d0*ah1(i)/27.d0)*dh1Zxx*dh1Zzz
     .       -(1.d0/3.d0)*dh1Zzz*dh2Zxx-(1.d0/3.d0)*dh1Zxx*dh2Zzz
      dat32(3,2,i)=dat32(3,2,i)+(12.d0*ah1(i)/27.d0)*dh1Zyy*dh1Zzz
     .       -(1.d0/3.d0)*dh1Zzz*dh2Zyy-(1.d0/3.d0)*dh1Zyy*dh2Zzz
      dat32(3,3,i)=dat32(3,3,i)+(12.d0*ah1(i)/27.d0)*dh1Zzz*dh1Zzz
     .       -(1.d0/3.d0)*dh1Zzz*dh2Zzz-(1.d0/3.d0)*dh1Zzz*dh2Zzz
      dat32(3,4,i)=dat32(3,4,i)+(12.d0*ah1(i)/27.d0)*dh1Zxy*dh1Zzz
     .       -(1.d0/3.d0)*dh1Zzz*dh2Zxy-(1.d0/3.d0)*dh1Zxy*dh2Zzz
      dat32(3,5,i)=dat32(3,5,i)+(12.d0*ah1(i)/27.d0)*dh1Zyz*dh1Zzz
     .       -(1.d0/3.d0)*dh1Zzz*dh2Zyz-(1.d0/3.d0)*dh1Zyz*dh2Zzz
      dat32(3,6,i)=dat32(3,6,i)+(12.d0*ah1(i)/27.d0)*dh1Zzx*dh1zzz
     .       -(1.d0/3.d0)*dh1Zzz*dh2Zzx-(1.d0/3.d0)*dh1Zzx*dh2Zzz
c
      dat32(4,1,i)=dat32(4,1,i)+(12.d0*ah1(i)/27.d0)*dh1Zxx*dh1Zxy
     .       -(1.d0/3.d0)*dh1Zxy*dh2Zxx-(1.d0/3.d0)*dh1Zxx*dh2Zxy
      dat32(4,2,i)=dat32(4,2,i)+(12.d0*ah1(i)/27.d0)*dh1Zyy*dh1Zxy
     .       -(1.d0/3.d0)*dh1Zxy*dh2Zyy-(1.d0/3.d0)*dh1Zyy*dh2Zxy
      dat32(4,3,i)=dat32(4,3,i)+(12.d0*ah1(i)/27.d0)*dh1Zzz*dh1Zxy
     .       -(1.d0/3.d0)*dh1Zxy*dh2Zzz-(1.d0/3.d0)*dh1Zzz*dh2Zxy
      dat32(4,4,i)=dat32(4,4,i)+(12.d0*ah1(i)/27.d0)*dh1Zxy*dh1Zxy
     .       -(1.d0/3.d0)*dh1Zxy*dh2Zxy-(1.d0/3.d0)*dh1Zxy*dh2Zxy
      dat32(4,5,i)=dat32(4,5,i)+(12.d0*ah1(i)/27.d0)*dh1Zyz*dh1Zxy
     .       -(1.d0/3.d0)*dh1Zxy*dh2Zyz-(1.d0/3.d0)*dh1Zyz*dh2Zxy
      dat32(4,6,i)=dat32(4,6,i)+(12.d0*ah1(i)/27.d0)*dh1Zzx*dh1zxy
     .       -(1.d0/3.d0)*dh1Zxy*dh2Zzx-(1.d0/3.d0)*dh1Zzx*dh2Zxy
c
      dat32(5,1,i)=dat32(5,1,i)+(12.d0*ah1(i)/27.d0)*dh1Zxx*dh1Zyz
     .       -(1.d0/3.d0)*dh1Zyz*dh2Zxx-(1.d0/3.d0)*dh1Zxx*dh2Zyz
      dat32(5,2,i)=dat32(5,2,i)+(12.d0*ah1(i)/27.d0)*dh1Zyy*dh1Zyz
     .       -(1.d0/3.d0)*dh1Zyz*dh2Zyy-(1.d0/3.d0)*dh1Zyy*dh2Zyz
      dat32(5,3,i)=dat32(5,3,i)+(12.d0*ah1(i)/27.d0)*dh1Zzz*dh1Zyz
     .       -(1.d0/3.d0)*dh1Zyz*dh2Zzz-(1.d0/3.d0)*dh1Zzz*dh2Zyz
      dat32(5,4,i)=dat32(5,4,i)+(12.d0*ah1(i)/27.d0)*dh1Zxy*dh1Zyz
     .       -(1.d0/3.d0)*dh1Zyz*dh2Zxy-(1.d0/3.d0)*dh1Zxy*dh2Zyz
      dat32(5,5,i)=dat32(5,5,i)+(12.d0*ah1(i)/27.d0)*dh1Zyz*dh1Zyz
     .       -(1.d0/3.d0)*dh1Zyz*dh2Zyz-(1.d0/3.d0)*dh1Zyz*dh2Zyz
      dat32(5,6,i)=dat32(5,6,i)+(12.d0*ah1(i)/27.d0)*dh1Zzx*dh1Zyz
     .       -(1.d0/3.d0)*dh1Zyz*dh2Zzx-(1.d0/3.d0)*dh1Zzx*dh2Zyz
c
      dat32(6,1,i)=dat32(6,1,i)+(12.d0*ah1(i)/27.d0)*dh1Zxx*dh1Zzx
     .       -(1.d0/3.d0)*dh1Zzx*dh2Zxx-(1.d0/3.d0)*dh1Zxx*dh2Zzx
      dat32(6,2,i)=dat32(6,2,i)+(12.d0*ah1(i)/27.d0)*dh1Zyy*dh1Zzx
     .       -(1.d0/3.d0)*dh1Zzx*dh2Zyy-(1.d0/3.d0)*dh1Zyy*dh2Zzx
      dat32(6,3,i)=dat32(6,3,i)+(12.d0*ah1(i)/27.d0)*dh1Zzz*dh1Zzx
     .       -(1.d0/3.d0)*dh1Zzx*dh2Zzz-(1.d0/3.d0)*dh1Zzz*dh2Zzx
      dat32(6,4,i)=dat32(6,4,i)+(12.d0*ah1(i)/27.d0)*dh1Zxy*dh1Zzx
     .       -(1.d0/3.d0)*dh1Zzx*dh2Zxy-(1.d0/3.d0)*dh1Zxy*dh2Zzx
      dat32(6,5,i)=dat32(6,5,i)+(12.d0*ah1(i)/27.d0)*dh1Zyz*dh1Zzx
     .       -(1.d0/3.d0)*dh1Zzx*dh2Zyz-(1.d0/3.d0)*dh1Zyz*dh2Zzx
      dat32(6,6,i)=dat32(6,6,i)+(12.d0*ah1(i)/27.d0)*dh1Zzx*dh1zzx
     .       -(1.d0/3.d0)*dh1Zzx*dh2Zzx-(1.d0/3.d0)*dh1Zzx*dh2Zzx
c
c     second derivatives of a
c
c      first term ( up to 'x')
c
      afac1=2.d0*dsqrt(at2(i)**3/27.d0)
      afac2=4.d0*dsqrt((at2(i)**3/27.d0)**3)
      afac1=max(afac1,1.e-16)
      afac2=max(afac2,1.e-16)
      afac1=ONE/afac1
      afac2=at3(i)/afac2
      ahelp=dsqrt(at2(i)**3/27.0d0)
      ahelp3=ONE-(at3(i)/TWO/ahelp)**2
      ahelp4=max(ahelp3,1.e-16)
      afac3=ahelp4**(-THREE_HALF)
      afac4=ahelp4**(-HALF)
c
c     similar zeroing out as in the first derivatives
c     unlike with the first derivatives, here it makes a difference
c
       if (ahelp3.lt.1.e-16) then
       afac3=ZERO
       afac4=ZERO
       endif
c
      !write (*,*) 'afac1, afac2 ', afac1, afac2
      do 106 k=1,6
      do 116 j=1,6
      da2(k,j,i)=(afac1*dat3(j,i)-afac2*(at2(i)**2/9.d0)*dat2(j,i))
     .          *(afac1*dat3(k,i)-afac2*(at2(i)**2/9.d0)*dat2(k,i)) 
           ! write (*,*) ' j, k , da2(j,k,i) ',  j, k , da2(j,k,i), dat3(j,i)
 116  continue
 106  continue
c
      afac33=at3(i)/2.d0/dsqrt(at2(i)**3/27.d0)
c     afac33=1 in uniaxial and afac33=-1 in equibiaxial
c
      do 117 k=1,6
      do 118 j=1,6
      da2(k,j,i)=-afac3*afac33*da2(k,j,i)

 118  continue
 117  continue
c
c     second term 'x'
c
      xcap=dsqrt(at2(i)**3/27.d0)
      xcap=max(xcap,1.e-16)
      xcap=1.d0/xcap
      dxcap=-at2(i)**2.d0*xcap**3/18.d0
      x1=xcap/2.d0
      x2=-xcap**3*at3(i)*at2(i)**2/36.d0
      x3=dxcap/2.d0
      x4=-xcap**2*dxcap*at3(i)*at2(i)**2/12.d0
      x5=-xcap**3*at2(i)**2/36.d0
      x6=-xcap**3*at3(i)*at2(i)/18.d0
c
      do 109 j=1,6
      do 110 k=1,6
      dx2(j,k,i)=x1*dat32(j,k,i)
     .          +x2*dat22(j,k,i)
     .          +x3*dat2(j,i)*dat3(k,i)
     .          +x4*dat2(j,i)*dat2(k,i)
     .          +x5*dat3(j,i)*dat2(k,i)
     .          +x6*dat2(j,i)*dat2(k,i)    
       !write (*,*) ' j, k , dat32(j,k,i) ',  j, k , dat32(j,k,i)       
      !write (*,*) ' j, k ,  dat3 (k,i)   ',  j, k , dat3(k,i)                 
 110  continue
 109  continue
c
      !print*, 'afac4 ' , afac4
      do 107 j=1,6
      do 108 k=1,6
      da2(j,k,i)=da2(j,k,i)-afac4*dx2(j,k,i)
 108  continue
 107  continue
c
c
c     second derivatives of principal strains
c     makes use of first derivatives of T2 : dat2
c     makes use of first derivatives of a : da
c     makes use of second derivatives of T2 : dat22
c     makes use of second derivatives of a : da2
c
      do 100 l=1,3
      ac=dcos(-aa4(i)/3.d0+(l-1)*2.d0*pi/3.d0)
      as=dsin(-aa4(i)/3.d0+(l-1)*2.d0*pi/3.d0) 
      do 101 j=1,6
      do 102 k=1,6
      t1=(-2.d0/12.d0/aa2(i)**3/3.d0)*ac*dat2(j,i)*dat2(k,i)
      t2=(2.d0/6.d0/aa2(i))          *ac*dat22(j,k,i)

      t3=(2.d0/18.d0/aa2(i))         *as*dat2(k,i)*da(j,i)
      t4=(2.d0/6.d0/aa2(i)/3.d0)     *as*dat2(j,i)*da(k,i)

      t5=(-2.d0/9.d0)*aa2(i)         *ac*da(j,i)*da(k,i)
      t6=(2.d0/3.d0) *aa2(i)         *as*da2(j,k,i)
      dsecond(l,j,k,i)=t1+t2+t3+t4+t5+t6
      !write (*,*) 't1, t2 , t3  ', t1, t2 , t3 
      !write (*,*) 't4, t5 , t6  ', t4, t5 , t6
      !write (*,*) ' j, k , da2(j,k,i) ',  j, k , da2(j,k,i)
 102  continue
 101  continue
 100  continue      
c
c
c
      enddo
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of calculation of
c     principal strain derivatives
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute static principal 2PK stresses
c     compute current tangent to stress-strain curve
c     in general : use the load curve with lowest strain rate
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   
c
      STRAINRATE(1:NEL,1)= ZERO
      STRAINRATE(1:NEL,2)= ZERO
      STRAINRATE(1:NEL,3)= ZERO
       !do i=1,NEL
      ! print*, 'JAC(I) ', i, JAC(I)
      ! enddo
      !do l=1,3
      
c
      do i=1,NEL
c
c    volumetric strain from Jacobian
c    positive in compression, for instance 0.8 means 80% compression
c    to be used as abcissa for the TABLE3D
c
      volstr(I) = ONE-JAC(I)
      !print*, 'JAC(I) ', i, JAC(I)
      !print*, 'volstr(I) ', i, volstr(I)
 
c     principal stretches from principal GL strains
c
       alambda(I,1)=dsqrt(2.d0*Z1(i)+1.d0)
       alambda(I,2)=dsqrt(2.d0*Z2(i)+1.d0)
       alambda(I,3)=dsqrt(2.d0*Z3(i)+1.d0)
c
c     principal engineering strain positive in compression
c     to be used as abcissa for the load curves
c
      EPSILON(I,1) = ONE - alambda(I,1)  
      EPSILON(I,2) = ONE - alambda(I,2)  
      EPSILON(I,3) = ONE - alambda(I,3)  
c
c     engineering strain rate, computed from true strain rate
c     and we take the absolute value ( so linear strain rate definition )
c     to be used as abcissa for the TABLE3D (SECOND ABCISSA)
c

      STRAINRATE(I,1)= abs(DRATE(i,1))
      STRAINRATE(I,2)= abs(DRATE(i,2))
      STRAINRATE(I,3)= abs(DRATE(i,3))

       !print*, 'engirate ' , 'STRAINRATE(I,1), STRAINRATE(I,2)',I, STRAINRATE(I,1), STRAINRATE(I,2)
       !print*, 'engirate ' , 'STRAINRATE(I,3)',STRAINRATE(I,3)
      !!! we only need true strain rate not the engineering one
c
      enddo !i=1,NEL
c
c     2 stress lookups are done now, static and dynamic
c     both use the TABLE2D that corresponds to the current value of volstr
c  

c     principal engineering stress positive in compression
c     this is from the quasistatic load curve ( lowest strain rate )
      XVEC(1:NEL,1)   = EPSILON(1:NEL,1)
      XVEC(1:NEL,2)   = ZERO           ! quasistatic
      XVEC(1:NEL,3)   = VOLSTR(1:NEL)

      IPOS(1:NEL,1) = 1
      IPOS(1:NEL,2) = 1
      IPOS(1:NEL,3) = 1            
      CALL TABLE_VINTERP(FUNC_SIG,NEL,IPOS,XVEC,FSIG1,DSIG1)

      XVEC(1:NEL,1)   = EPSILON(1:NEL,2)
      XVEC(1:NEL,2)   = ZERO           ! quasistatic
      XVEC(1:NEL,3)   = VOLSTR(1:NEL)

      IPOS(1:NEL,1) = 1
      IPOS(1:NEL,2) = 1
      IPOS(1:NEL,3) = 1            
      CALL TABLE_VINTERP(FUNC_SIG,NEL,IPOS,XVEC,FSIG2,DSIG2)

      XVEC(1:NEL,1)   = EPSILON(1:NEL,3)
      XVEC(1:NEL,2)   = ZERO           ! quasistatic
      XVEC(1:NEL,3)   = VOLSTR(1:NEL)

      IPOS(1:NEL,1) = 1
      IPOS(1:NEL,2) = 1
      IPOS(1:NEL,3) = 1            
      CALL TABLE_VINTERP(FUNC_SIG,NEL,IPOS,XVEC,FSIG3,DSIG3)

      do i=1,NEL
c
        sigma(I,1) = FSIG1(I) * SCAL
        slope(I,1) = DSIG1(I) * SCAL

        sigma(I,2) = FSIG2(I) * SCAL
        slope(I,2) = DSIG2(I) * SCAL

        sigma(I,3) = FSIG3(I) * SCAL
        slope(I,3) = DSIG3(I) * SCAL

      enddo !i=1,NEL
c
c     principal engineering stress positive in compression
c     this is from the complete rate dependent TABLE2D
c     so 'engirate' has to come in as an additional argument
c     DIRECTION 1
      XVEC(1:NEL,1)   = EPSILON(1:NEL,1)
      XVEC(1:NEL,2)   = STRAINRATE(1:NEL,1) !  *XFAC
      XVEC(1:NEL,3)   = VOLSTR(1:NEL)
!       
      IPOS(1:NEL,1) = 1
      IPOS(1:NEL,2) = 1
      IPOS(1:NEL,3) = 1            
      CALL TABLE_VINTERP(FUNC_SIG,NEL,IPOS,XVEC,FSIG1,DSIG1)

c     DIRECTION 2

      XVEC(1:NEL,1)   = EPSILON(1:NEL,2)
      XVEC(1:NEL,2)   = STRAINRATE(1:NEL,2) !  *XFAC
      XVEC(1:NEL,3)   = VOLSTR(1:NEL)
!       
      IPOS(1:NEL,1) = 1
      IPOS(1:NEL,2) = 1
      IPOS(1:NEL,3) = 1            
      CALL TABLE_VINTERP(FUNC_SIG,NEL,IPOS,XVEC,FSIG2,DSIG2)
      
c     DIRECTION 3

      XVEC(1:NEL,1)   = EPSILON(1:NEL,3)
      XVEC(1:NEL,2)   = STRAINRATE(1:NEL,3) !  *XFAC
      XVEC(1:NEL,3)   = VOLSTR(1:NEL)
!       
      IPOS(1:NEL,1) = 1
      IPOS(1:NEL,2) = 1
      IPOS(1:NEL,3) = 1            
      CALL TABLE_VINTERP(FUNC_SIG,NEL,IPOS,XVEC,FSIG3,DSIG3)

      do i=1,NEL

      dynsigma(I,1) = FSIG1(I) *  SCAL
      dynsigma(I,2) = FSIG2(I) *  SCAL
      dynsigma(I,3) = FSIG3(I) *  SCAL
c
c     principal 2PK stress positive in tension
c
      sik(i,1)=-sigma(I,1)/alambda(I,1)
      sik(i,2)=-sigma(I,2)/alambda(I,3)
      sik(i,3)=-sigma(I,3)/alambda(I,3)
c
c     incremental principal stiffnesses
c
      dik(i,1) = slope(i,1)/alambda(I,1)**2+sigma(I,1)/alambda(I,1)**3
      dik(i,2) = slope(i,2)/alambda(I,2)**2+sigma(I,2)/alambda(I,2)**3
      dik(i,3) = slope(i,3)/alambda(I,3)**2+sigma(I,3)/alambda(I,3)**3

c
c     viscous overstress in terms of 2PK
c     positive in tension
c
      drate(i,1)=-(dynsigma(I,1)-sigma(I,1))/alambda(I,1)
      drate(i,2)=-(dynsigma(I,2)-sigma(I,2))/alambda(I,2)
      drate(i,3)=-(dynsigma(I,3)-sigma(I,3))/alambda(I,3)

      !write(*,*)' sik(i,l) , dik(i,l) ', sik(i,l)  , dik(i,l)


      enddo
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute incremental stiffness matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      do l=1,3
      do k=1,6
      do j=1,6
      do i=1,NEL
      cijkl(i,j,k)=ZERO
      enddo
      enddo
      enddo
      enddo
c
      do k=1,6
      do j=1,6
      do i=1,NEL
      do l=1,3
      cijkl(i,j,k)=cijkl(i,j,k)
     .  + dfirst(l,j,i)*dfirst(l,k,i)*dik(i,l)
     .  + dsecond(l,j,k,i)*sik(i,l)
       !write(*,*)' dfirst(l,j,i),dfirst(l,k,i) ', dfirst(l,j,i),dfirst(l,k,i) 
       !write(*,*)'   j,k,l ,      dsecond(l,j,k,i) ', j,k,l ,dsecond(l,j,k,i)
   
      enddo
      enddo
      enddo
      enddo
c
c     factor 1/4 for Voigt notation with gammas ( 2 gammas)
c     factor 1/2 for Voigt notation with gammas ( 1 gamma)
c
      do k=4,6
      do j=4,6
      do i=1,NEL
      cijkl(i,j,k)=cijkl(i,j,k)/4.
      enddo
      enddo
      enddo
c
      do k=1,3
      do j=4,6
      do i=1,NEL
      cijkl(i,j,k)=cijkl(i,j,k)/2.
      enddo
      enddo
      enddo
      do k=4,6
      do j=1,3
      do i=1,NEL
      cijkl(i,j,k)=cijkl(i,j,k)/2.
      enddo
      enddo
      enddo
c
c
c
      return
      end
