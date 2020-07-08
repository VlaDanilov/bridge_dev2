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
 
      function A51mpppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A51mpppp
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5
C----Eq.(4) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(-i)*(-1)
C--- to give (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      A51mpppp=-((s(j2,j3)+s(j3,j4)+s(j4,j5))*zb(j2,j5)**2
     & -zb(j2,j4)*za(j4,j3)*zb(j3,j5)*zb(j2,j5)
     & -zb(j1,j2)*zb(j1,j5)/(za(j1,j2)*za(j1,j5))
     & *(za(j1,j2)**2*za(j1,j3)**2*zb(j2,j3)/za(j2,j3)
     &  +za(j1,j3)**2*za(j1,j4)**2*zb(j3,j4)/za(j3,j4)
     &  +za(j1,j4)**2*za(j1,j5)**2*zb(j4,j5)/za(j4,j5)))
     & /(3._dp*zb(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))
      return
      end 
