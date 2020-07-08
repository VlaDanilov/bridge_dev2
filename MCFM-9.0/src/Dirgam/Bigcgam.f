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
 
      function BigCgam(i1,i2,i3,i4,i5,i)
      implicit none
      include 'types.f'
      real(dp):: BigCgam
      
CCCCCC Matrix element squared for
C     qbar(-p1)+q(-p2)=g(p3)+g(p4)+gamma(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer:: i1,i2,i3,i4,i5,i
      BigCgam=
     & (xn*(s(i1,i3)*s(i2,i4)+s(i2,i3)*s(i1,i4))/s(i4,i3)-s(i1,i2)/xn)
     .*((s(i2,i5)**2+s(i1,i5)**2)/(s(i1,i4)*s(i2,i4)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i4)**2+s(i4,i1)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i3)**2+s(i1,i3)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i4)*s(i2,i4)))

      BigCgam=32._dp*esq*gsq**2*Q(i)**2*BigCgam

      return
      end
