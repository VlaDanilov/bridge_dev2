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
 
      function aqqb_wbb(i1,i2,i5,i6,i3,i4)
      implicit none
      include 'types.f'
      complex(dp):: aqqb_wbb
      
C---This is the amplitude for 
C---q_L(p1)+q_L(p6) --> q_L(p2)+q_L(p5)+W(l(p3)+antilepton(p4))
c---with no couplings included
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
C---'prods.f' includes both the declaration and the common for s and za,zb
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: t2
      real(dp):: s234,s256,prop
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
      s256=s(i2,i6)+s(i2,i5)+s(i5,i6)
      prop=s(i5,i6)*sqrt((s(i3,i4)-wmass**2)**2+(wmass*wwidth)**2)
      aqqb_wbb=
     & +za(i3,i2)*zb(i6,i1)*t2(i4,i2,i3,i5)/(prop*s234)
     & +za(i5,i2)*zb(i4,i1)*t2(i6,i2,i5,i3)/(prop*s256)
      return 
      end
