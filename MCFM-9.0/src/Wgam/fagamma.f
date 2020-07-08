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
 
      function fagamma(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fagamma
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: L0,L1,Lsm1
      fagamma=
     & +za(i1,i3)**2/(za(i1,i5)*za(i3,i4)*za(i2,i5))
     & *Lsm1(-s(i2,i5),-s(i3,i4),-s(i1,i2),-s(i3,i4))
     & +0.5_dp*za(i1,i2)*za(i3,i5)*zb(i4,i5)*zb(i2,i5)
     & /(za(i2,i5)*s(i1,i2)**2)*L1(-s(i3,i4),-s(i1,i2))
     & -1.5_dp*za(i1,i3)*zb(i4,i5)/(za(i2,i5)*s(i1,i2))
     & *L0(-s(i3,i4),-s(i1,i2))
     & -0.5_dp*zb(i2,i4)*zb(i4,i5)/(zb(i1,i2)*zb(i3,i4)*za(i2,i5))
      return
      end
