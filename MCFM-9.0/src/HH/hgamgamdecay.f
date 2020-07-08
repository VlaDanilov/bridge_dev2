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
 
      subroutine hgamgamdecay(p,j,k,hdecay)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'couple.f'
      include 'kpart.f'
      include 'first.f'
      integer:: j,k
      complex(dp):: Iw,Iq,Ftriangle
      real(dp):: p(mxpart,4),prefac,mh,hdecay
      real(dp):: x_t,x_b,x_w,x,mt_eff,mb_eff,massfrun
      save mt_eff,mb_eff
!$threadprivate(mt_eff,mb_eff)
C---statement functions
      Iq(x)=cplx1(4._dp*x)*(ctwo+cplx1(4._dp*x-1._dp)*Ftriangle(x))
      Iw(x)=-ctwo*(cplx1(6._dp*x+1._dp)
     & +cplx1(6._dp*x*(2._dp*x-1._dp))*Ftriangle(x))
C---end statement functions
      mh=sqrt((p(j,4)+p(k,4))**2
     &       -(p(j,1)+p(k,1))**2
     &       -(p(j,2)+p(k,2))**2
     &       -(p(j,3)+p(k,3))**2)

      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mb_eff=massfrun(mb_msbar,hmass,amz,1)
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mb_eff=massfrun(mb_msbar,hmass,amz,2)
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

C---Total width multiplied by a factor of 16*pi*mh 
C---to get matrix element squared.
c---maybe it would be better to add esq at a higher scale.
      prefac=(esq/(4._dp*pi))**2*Gf*mh**4/(8._dp*rt2*pi**2)
      x_b=(mb_eff/mh)**2
      x_t=(mt_eff/mh)**2
      x_w=(wmass/mh)**2
      hdecay=prefac
     & *abs(xn*(Q(1)**2*Iq(x_b)+Q(2)**2*Iq(x_t))+Iw(x_w))**2
      return
      end
