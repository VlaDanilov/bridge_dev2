!  Copyright (C) 2019, respective authors of MCFM.
!
!  This program is free software: you can redistribute it and/or modify it under
!  the terms of the GNU General Public License as published by the Free Software
!  Foundation, either version 3 of the License, or (at your option) any later
!  version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!  PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program. If not, see <http://www.gnu.org/licenses/>
 
      subroutine qqb_tautau(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: JMC                  
*     May, 1999.                           
*     calculate the element squared 
*     for the process in terms of p     
*     q(-p1) +qbar(-p2) -->  tau^+(nubar_tau(p6)+nu_e(p7)+e+(p8))               
*                           +tau^-(nu_tau(p5)+e-(p3)+nubar_e(p4))             
*     This routine is nothing more than a wrapper for 
*     Kleiss and Stirling
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'first.f'
      include 'mpicommon.f'
      integer:: j,nu
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),pks(4,10)
      real(dp):: RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP
      real(dp):: GW_KS,GS_KS,fac
      complex(dp):: qqbpp,qqbpm,qqbmp,qqbmm
      complex(dp):: qbqpp,qbqpm,qbqmp,qbqmm,prop12
      real(dp):: s12,dotks
      COMMON/COUPS/GW_KS,GS_KS
      COMMON/PARS/RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP
      COMMON/MOM/PKS
!$omp threadprivate(/COUPS/,/PARS/,/MOM/)  

c---- Fill common blocks for Kleiss and Stirling   
      if (first) then
         first=.false.
         rmw=wmass
         rgw=wwidth
         rmt=mtau
         rgt=tauwidth
         rmb=zip
         rmtlo=zip
         rmtup=zip
         GS_KS=sqrt(gsq)
         GW_KS=sqrt(gwsq/eight)
!$omp master
         if (rank.eq.0) then
         write(6,*) 
         write(6,*) 'rmw',rmw
         write(6,*) 'rmt',rmt
         write(6,*) 'rgw',rgw
         write(6,*) 'rgt',rgt
         write(6,*) 'GS_KS',GS_KS
         write(6,*) 'GW_KS',GW_KS
         write(6,*) 
         endif
!$omp end master
      endif

c---step one change momentum notation into the notation of Kleiss and Stirling
c----My notation
*     q(-p1) +qbar(-p2) -->  tau^+(nubar_tau(p6)+nu_e(p7)+e+(p8))
*                           +tau^-(nu_tau(p5)+e-(p3)+nubar_e(p4))
c----KS notation
c---qbar(p1)+q(p2) = tau^+(nubar_tau(p3)+e^-(p4)+nubar_e(p5))
c                   +tau^-(nu_tau(p6)+nu_e(p7)+e^+(p8))

      do nu=1,4
      pks(nu,1)=-p(2,nu)
      pks(nu,2)=-p(1,nu)
      pks(nu,3)=p(6,nu)
      pks(nu,4)=p(7,nu)
      pks(nu,5)=p(8,nu)
      pks(nu,6)=p(5,nu)
      pks(nu,7)=p(3,nu)
      pks(nu,8)=p(4,nu)
      enddo


C-----Matrix elements of Kleiss and Stirling include all averaging

      call tautauww(1,2,qqbpp,qqbpm,qqbmp,qqbmm,fac)
      qbqpp=qqbmp
      qbqpm=qqbmm
      qbqmp=qqbpp
      qbqmm=qqbpm
C----set all elements to zero
      msq(:,:)=zip
      
      s12=two*dotks(1,2)
      prop12=s12/cplx2((s12-zmass**2),zmass*zwidth)
      fac=fac*esq**2
C---fill qb-q and q-qb elements
      do j=-nf,nf
      if (j < 0) then
        msq(j,-j)=fac*(abs(qbqpp*(R(-j)*re*prop12-Q(j))
     &                    +qbqpm*(R(-j)*le*prop12-Q(j)))**2
     &                +abs(qbqmp*(L(-j)*re*prop12-Q(j))
     &                    +qbqmm*(L(-j)*le*prop12-Q(j)))**2)
      elseif (j > 0) then
        msq(j,-j)=fac*(abs(qqbpp*(R(+j)*re*prop12-Q(j))
     &                    +qqbpm*(R(+j)*le*prop12-Q(j)))**2
     &                +abs(qqbmp*(L(+j)*re*prop12-Q(j))
     &                    +qqbmm*(L(+j)*le*prop12-Q(j)))**2)
      endif
      enddo
      
      return
      end
