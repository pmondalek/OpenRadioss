! copyright>        openradioss
! copyright>        copyright (c) 1986-2024 altair engineering inc.
! copyright>
! copyright>        this program is free software: you can redistribute it and/or modify
! copyright>        it under the terms of the gnu affero general public license as published by
! copyright>        the free software foundation, either version 3 of the license, or
! copyright>        (at your option) any later version.
! copyright>
! copyright>        this program is distributed in the hope that it will be useful,
! copyright>        but without any warranty; without even the implied warranty of
! copyright>        merchantability or fitness for a particular purpose.  see the
! copyright>        gnu affero general public license for more details.
! copyright>
! copyright>        you should have received a copy of the gnu affero general public license
! copyright>        along with this program.  if not, see <https://www.gnu.org/licenses/>.
! copyright>
! copyright>
! copyright>        commercial alternative: altair radioss software
! copyright>
! copyright>        as an alternative to this open-source version, altair also offers altair radioss
! copyright>        software under a commercial license.  contact altair to discuss further if the
! copyright>        commercial version may interest you: https://www.altair.com/radioss/.
! chd|====================================================================
! chd|  sigeps169_connect            source/materials/mat/mat169/sigeps169_connect.f
! chd|-- called by -----------
! chd|        suser43                       source/elements/solid/sconnect/suser43.f
! chd|-- calls ---------------
! chd|        finter                        source/tools/curve/finter.f   
! chd|====================================================================
! ======================================================================================================================
      module sigeps169_connect_mod
      contains
! ======================================================================================================================
! material for cohesive element, elastic in normal direction, elastoplastic in shear, with coupled damage
! ======================================================================================================================         
        subroutine sigeps169_connect(                                          &
          nel     ,time, uparam,                                               &
          nuparam ,stifm   ,                                                   &
          area    , off     ,nuvar  ,uvar  ,ipg      ,                         &
          depszz  ,depsyz   ,depszx ,epszz   ,epsyz    ,epszx   ,              &                             
          sigozz  ,sigoyz  ,sigozx  ,signzz  ,signyz ,signzx  ,                &
          pla     ,iout    ,jsms    ,dmg     ,ngl   ,dmels ,                   &
          idtmins,dtfacs   ,dtmins  , thick0)                                        

!-----------------------------------------------
!   m o d u l e s
!-----------------------------------------------
      use constant_mod , only : one , zero,two,em20,four,five,half
!-----------------------------------------------
      implicit none
! ----------------------------------------------
#include "my_real.inc"
!-----------------------------------------------
!   d u m m y   a r g u m e n t s
!-----------------------------------------------       
      integer ,intent(in)  :: nel,nuparam,nuvar,jsms,ipg,iout

      integer ,dimension(nel)       ,intent(in)    :: ngl

      integer ,intent(in) :: idtmins
      my_real ,intent(in) :: dtfacs
      my_real ,intent(in) :: dtmins

      my_real ,intent(in)  :: time
      my_real ,dimension(nel)  ,intent(inout) :: off,area,pla,dmels
      my_real ,dimension(nel)  ,intent(in)    :: depszz,depsyz,depszx,epszz,epsyz,epszx,thick0
      my_real ,dimension(nel)  ,intent(in)    :: sigozz  ,sigoyz  ,sigozx  
      my_real ,dimension(nel)  ,intent(out)   :: signzz,signyz,signzx
      my_real ,dimension(nel)  ,intent(inout) :: dmg,stifm
      my_real ,dimension(nuparam)   ,intent(in)    :: uparam
      my_real ,dimension(nel,nuvar) ,intent(inout) :: uvar

!-----------------------------------------------
!   l o c a l   v a r i a b l e s
!-----------------------------------------------
        integer      i,ii,k,iterk,nindf,tempi,nindxdn,nindxdsh
        integer  iel
        integer ,dimension(nel)   ::      indxdsh,indxdn,      indf

        my_real        young,nu, wave,tenmax,gcten,shrmax, gcshr,shrp, sht_sl,pwrt,pwrs ,taumax
        my_real        dlam, dpla_dlam,shear,an,as,dp,g1,g2,sigeq,dtb,norm,nxz,nyz,d0fn,dfn,d0fs,dfs
        my_real, dimension(nel) :: dstrnz, dstrnxz, dstrnyz,strs_tr_sh,fyld,strn_pl,dstr_sh,epsn0
        my_real, dimension(nel) :: normzz, normxz, normyz, df,dpzx,dpyz,dpzz,ddpla,dpla,dmg_n,dmg_s
        my_real, dimension(nel) :: fdam_n,fdam_s,thick,stf,tempr,normyz_norm,normxz_norm,normzz_norm,eps_sh0,eps_sh

!--------------------------------------------------------
     


      young    =   uparam(1) 
      shear    =   uparam(2)
      nu       =   uparam(3) 
      tenmax   =   uparam(4) 
      gcten    =   uparam(5) 
      shrmax   =   uparam(6) 
      gcshr    =   uparam(7) 
      shrp     =   uparam(8) 
      sht_sl   =   uparam(9) 
      dfn      =   uparam(10)
      dfs      =   uparam(11)
      dp       =   uparam(12)
      pwrt     =   uparam(13)
      pwrs     =   uparam(14)
      iterk    = 3


       
      do i=1,nel 
        epsn0(i)   = uvar(i,3)        
        eps_sh0(i) = uvar(i,4)    
      enddo

      wave = young * (one-nu)/(one+nu)/(one-two*nu)
      stf(1:nel)     = young *  area(1:nel)                                       
      stifm(1:nel)   = stifm(1:nel)  + stf(1:nel)*off(1:nel)                                              
      ! omega = sqrt(2k/2*dmels), dt=2/omega, 2*dmels=dt**2 * 2k / 4
      if (idtmins==2 .and. jsms/=0) then
        dtb = (dtmins/dtfacs)**2
        do iel=1,nel                                                 
          dmels(iel)=max(dmels(iel),half*dtb*stf(iel)*off(iel))
        enddo                                                        
      end if                                                        

      an      = tenmax/gcten/two
      as      = shrmax/two/gcshr

      

        
      do i=1,nel  
        dmg_n(i) = uvar(i,1)
        dmg_s(i) = uvar(i,2)        
        fdam_n(i) = one - dmg_n(i)
        fdam_s(i) = one - dmg_s(i)
        wave = young * (one-nu)/(one+nu)/(one-two*nu)
      ! write(*,*) ' dmg_n       = ',i,  dmg_n(i) 
      enddo                                                        


      dpla(1:nel) = zero
     
      do i = 1,nel
        if (off(i) < 0.001)  off(i) = zero
        if (off(i) < one)    off(i) = off(i)*four/five
        if (off(i) == one) then
         
          signzz(i)     = sigozz(i)/fdam_n(i) + (depszz(i)  )*wave/thick0(i) !young is per unit length
          signzz(i)     = min(signzz(i),tenmax)
               ! write(*,*) ' signzz(i)  = ', i, signzz(i)

          signyz(i)     = sigoyz(i)/ fdam_s(i)+ (depsyz(i)  ) *shear/thick0(i)
          !signyz(i)     = min(signyz(i),shrmax)
          !signyz(i)     = min(signyz(i),sigoyz(i))
          signzx(i)     = sigozx(i)/ fdam_s(i)+ (depszx(i)  )*shear/thick0(i)
          !signzx(i)     = min(signzx(i),shrmax)
          !signzx(i)     = min(signzx(i),sigozx(i))
          strs_tr_sh(i) = sqrt( signyz(i)**2 + signzx(i)**2)
        endif
      enddo    
      !-----------------------------------------------
      !   compute yield function 
      !-----------------------------------------------
      do i=1,nel
            fyld(i) = (max(signzz(i),zero)/ tenmax)**pwrt +                         &
                      (strs_tr_sh(i) /(shrmax -sht_sl * signzz(i) ))**pwrs 
            fyld(i) =  fyld(i) - one        
            pla(i) = pla(i)
      enddo
         !-----------------------------------------------
         !   compute plasticity
         !-----------------------------------------------
      do i=1,nel
        !if ((g2**pwrs-one) > zero.and. off(i) == one )then    
        if (fyld(i) > zero.and. off(i) == one )then  
          g1       = max(signzz(i),zero)/ tenmax
          g2       = strs_tr_sh(i) /(shrmax - sht_sl * signzz(i))
          do k=1,iterk
             !computation of the normal to the yield surface  
            normxz(i) =  pwrs *  (g2 **(pwrs-1) )                           &
                          *  signzx(i)/max(strs_tr_sh(i) ,em20)             &
                         /(shrmax - sht_sl * signzz(i)) 

            normyz(i) =  pwrs * (g2 **(pwrs-1)   )                          &
                          * signyz(i)/max(strs_tr_sh(i),em20)               &
                          /(shrmax - sht_sl * signzz(i)) 
            !normzz(i) = zero ! pwrt  * g1 **(pwrt -one) / tenmax                  &
                   !        +pwrs * g2 **(pwrs-one)                          &
                   ! * (sht_sl* strs_tr_sh(i) / (shrmax - sht_sl * signzz(i))**2 )
            norm =  sqrt(normxz(i) **2 +  normyz(i)**2  )  !  + normzz(i)**2       
            normxz_norm(i) =   normxz(i)/norm
            normyz_norm(i) =   normyz(i)/norm       
            !normzz_norm(i) =   normzz(i)/norm
            !df/dlam = df/dsig * dsig/dlam = nij * celas* nij
            df(i) = -((shear/thick0(i)) *(normxz(i) **2 + normyz(i)**2 ))/  norm   
            df(i) = sign(max(abs(df(i)),em20),df(i)) 

            ! computation of the plastic multiplier
            dlam = -fyld(i)/df(i)
            ! plastic strains tensor increment  
            dpzx(i) = dlam * normxz_norm(i)
            dpyz(i) = dlam * normyz_norm(i)

            !dpzz(i) = dlam * normzz_norm(i)
            ! elasto-plastic stresses update   
            signyz(i)     =  signyz(i) - shear/thick0(i) * dpyz(i)
            signzx(i)     =  signzx(i) - shear/thick0(i) * dpzx(i)
            !signzz(i)     =  signzz(i) - wave  * dpzz(i)
            strs_tr_sh(i) =  sqrt( signyz(i)**2 + signzx(i)**2)
            g2        =  strs_tr_sh(i) /(shrmax - sht_sl * signzz(i) )
            !g1        = max(signzz(i),zero)/ tenmax
            sigeq     =  g1**pwrt + g2**pwrs
            !fyld(i)   =  g2**pwrs- one
            fyld(i)   =  sigeq- one   

            ddpla(i) = (signyz(i)*dpyz(i) + signzx(i)*dpzx(i) )/strs_tr_sh(i) !(shrmax - sht_sl * signzz(i))                
            dpla(i)  = max(zero, dpla(i) + ddpla(i) )  
            pla(i)   = pla(i)  +  ddpla(i)          
          enddo
        endif
        !print*, 'i, dpzx, dpyz ', i, dpzx(i),dpyz(i)
      enddo
   
      !  !  ! check if damage reached in normal direction
      nindxdn = 0
      nindf = 0
      do i=1,nel  
        if (signzz(i) >= tenmax .or. (fyld(i)>zero.and. signzz(i)>zero))  then 
        !if (signzz(i) >= tenmax .or. (dpla(i)>zero.and. signzz(i)>zero))  then         
          !compute damage (linear)   
          !dmg_n(i) =   dmg_n(i) + an *  depszz(i)   
          if(uvar(i,3) == zero) then
            epsn0(i)  = epszz(i)
            uvar(i,3) = epsn0(i)
            nindxdn = nindxdn+1
            indxdn(nindxdn) = i
          endif  
  
          dmg_n(i) =   (epszz(i) -epsn0(i)) /(dfn -epsn0(i))                   
          dmg_n(i) =   max(uvar(i,1),dmg_n(i)  )          

          if (dmg_n(i)>=one)then 
            dmg_n(i) = one
            off(i) = four/five
            nindf   = nindf + 1  
            indf(nindf) = i        
          endif 
        endif 
      enddo

      ! check si strs_tr_sh>shrmax et si dpla depasse le shrp (input) si oui calcul endommagement shear  
      nindxdsh = 0

      do i=1,nel
        strs_tr_sh(i) = sqrt(signyz(i)**2 + signzx(i)**2)
        dstr_sh(i)    = sqrt(depsyz(i)**2 + depszx(i)**2) 
        eps_sh(i)     = sqrt(epsyz(i)**2  + epszx(i)**2) 
        taumax        = shrmax - sht_sl * signzz(i)
        if (fyld(i) >zero) write(*,*) 'i, strs_tr_sh ', ngl(i) ,  strs_tr_sh(i)
        !if ((strs_tr_sh(i) >= taumax .or. dpla(i)>zero ).and. pla(i) >= dp )  then 
        if ((strs_tr_sh(i) >= taumax .or. fyld(i)>zero) .and. (eps_sh(i) - eps_sh0(i)) >= dp)then !.and. pla(i) >= dp 
          if(uvar(i,4)==zero) then 
            eps_sh0(i) = eps_sh(i)
            uvar(i,4)  = eps_sh0(i)
            nindxdsh = nindxdsh+1
            indxdsh(nindxdsh) = i
          endif  
          dmg_s(i) = (eps_sh(i) - eps_sh0(i)) / (dfs - eps_sh0(i))  
          !print*, '  eps_sh(i) , eps_sh0(i) ', eps_sh(i) , eps_sh0(i)
          !dmg_s(i) = min(one,max(zero,(dmg_s(i) + as * dstr_sh(i)  ) ) )

          if (dmg_s(i)>=one)then 
            dmg_s(i) = one
            if (dmg_n(i) /=one)then 
              off(i) = four/five
              nindf   = nindf + 1  
              indf(nindf) = i    
            endif  
          endif 
        endif      
      enddo
      

      do i=1,nel
        uvar(i,1) = dmg_n(i)  
        uvar(i,2) = dmg_s(i)  
        dmg(i) =max(dmg_n(i),dmg_s(i))

      enddo

      do i=1,nel
        strs_tr_sh(i) =  sqrt( signyz(i)**2 + signzx(i)**2)
        g2        =  strs_tr_sh(i) /(shrmax - sht_sl * signzz(i) )
        g1        = max(signzz(i),zero)/ tenmax
        sigeq     =  g1**pwrt + g2**pwrs
        uvar(i,5) = sigeq
                !write(*,*) ' signzz(i) dmgn = ', i, signzz(i), dmg_n(i)
        signzz(i)     =  signzz(i) * (one -  dmg_n(i))
        signyz(i)     =  signyz(i) * (one -  dmg_s(i))
        signzx(i)     =  signzx(i) * (one -  dmg_s(i))       
      enddo
   
!-------------------------
       if (nindxdn > 0) then
         do ii=1,nindxdn
           write(iout, 1000) ngl(indxdn(ii)),time
         enddo
       end if
       if (nindxdsh > 0) then
         do ii=1,nindxdsh
           write(iout, 1002) ngl(indxdsh(ii)),time
         enddo
       end if
       if (nindf > 0) then
         do ii=1,nindf
           write(iout, 1003) ngl(indf(ii)),ipg,time
           write(iout, 1004) epszz(indf(ii)), eps_sh(indf(ii))
         enddo
       end if        
!-----------------------------------------------------------------
    1000 format(1x,'START DAMAGE IN NORMAL DIRECTION IN CONNECTION ELEMENT NUMBER ',i10,1x,' AT TIME :',g11.4) 
    1002 format(1x,'START DAMAGE IN SHEAR DIRECTION IN CONNECTION ELEMENT NUMBER ',i10,1x,' AT TIME :',g11.4)     
    1003 format(1x,'FAILURE IN CONNECTION ELEMENT NUMBER ',i10,1x,'INTEGRATION POINT',i10,1x, 'AT TIME :',g11.4)     
    1004 format(1x,'ELONGATION IN NORMAL DIRECTION AT FAILURE ',g11.4,1x,'ELONGATION IN TANGENTIAL DIRECTION AT FAILURE',g11.4) 
!
!-----------------------------------------------------------------
     end subroutine sigeps169_connect
end module sigeps169_connect_mod             
