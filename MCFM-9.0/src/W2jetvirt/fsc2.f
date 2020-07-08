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
 
      function Fsc2(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Fsc2
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L1
      real(dp):: t  
      Fsc2= 
     .za(j3,j5)**2/(2._dp*za(j1,j2)*za(j2,j3)*za(j5,j6)*zb(j3,j4))+
     .zb(j2,j6)**2/(2._dp*za(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j5,j6))-
     .(za(j1,j3)*((L1(-t(j1,j2,j3),-s(j2,j3))*za(j2,j3)**2*zb(j1,j2)**2*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2)/s(j2,j3)**2+
     .L1(-s(j5,j6),-t(j1,j2,j3))*za(j3,j4)**2*zb(j4,j6)**2))/
     .(2._dp*za(j1,j2)*za(j2,j3)
     .*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))
     .*zb(j5,j6)*t(j1,j2,j3))-(zb(j2,j4)*
     .(L1(-s(j5,j6),-t(j2,j3,j4))*za(j1,j5)**2*zb(j1,j2)**2+
     .(L1(-t(j2,j3,j4),-s(j2,j3))*za(j3,j4)**2*zb(j2,j3)**2*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))**2)/s(j2,j3)**2))/
     .(2._dp*za(j5,j6)*zb(j2,j3)*zb(j3,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*t(j2,j3,j4))
      return
      end
