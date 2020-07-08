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
 
      function Bigbgam(i1,i2,i3,i4,i5,i)
      implicit none
      include 'types.f'
      real(dp):: Bigbgam
      
CC    Matrix element for q(-p1)+q(-p2)->q(p3)+q(p4)+gamma(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer:: i1,i2,i3,i4,i5,i
      real(dp):: eik,ss,ssp,tt,ttp,uu,uup,e12,e34,e13,e14,e23,e24
      e12=s(i1,i2)/(s(i1,i5)*s(i2,i5))
      e34=s(i3,i4)/(s(i3,i5)*s(i4,i5))
      e14=s(i1,i4)/(s(i1,i5)*s(i4,i5))
      e24=s(i2,i4)/(s(i2,i5)*s(i4,i5))
      e13=s(i1,i3)/(s(i1,i5)*s(i3,i5))
      e23=s(i2,i3)/(s(i2,i5)*s(i3,i5))
      eik=-e12+e23+e14-e34+e13+e24
      ss=s(i1,i2)
      ssp=s(i3,i4)
      tt=s(i1,i3)
      ttp=s(i2,i4)
      uu=s(i1,i4)
      uup=s(i2,i3)
      Bigbgam=8._dp*xn*CF*gsq**2*esq*Q(i)**2*eik
      Bigbgam=Bigbgam*
     & ((ss**2+ssp**2+uu**2+uup**2)/(tt*ttp)
     & +(ss**2+ssp**2+tt**2+ttp**2)/(uu*uup)
     & -(ss**2+ssp**2)*(ss*ssp-tt*ttp-uu*uup)/(xn*tt*ttp*uu*uup))

      return
      end
