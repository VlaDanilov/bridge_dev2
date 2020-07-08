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
 
      function A51ppppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      complex(dp):: A51ppppp
      integer:: j1,j2,j3,j4,j5
      complex(dp):: asym,sym
      asym=+zb(j1,j2)*za(j2,j3)*zb(j3,j4)*za(j4,j1)
     &     -za(j1,j2)*zb(j2,j3)*za(j3,j4)*zb(j4,j1)
      sym=cplx1(s(j1,j2)*s(j2,j3)+s(j2,j3)*s(j3,j4)
     &          +s(j3,j4)*s(j4,j5)+s(j4,j5)*s(j5,j1)+s(j5,j1)*s(j1,j2))

C----Eq.(4) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(-i)*(-1)
C--- to give (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      A51ppppp=-(sym+asym)
     % /(6._dp*za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1))
      return
      end 
