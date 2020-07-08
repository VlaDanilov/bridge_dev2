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
 
      function Fcc_qpqmgpgp(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Fcc_qpqmgpgp
      
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: L0,Lsm1,Lsm1_2me
      real(dp):: t  

      Fcc_qpqmgpgp=
     .((Lsm1(-s(j1,j4),-t(j1,j2,j4),-s(j1,j2),-t(j1,j2,j4))+
     .Lsm1_2me(t(j1,j3,j4),t(j1,j2,j4),s(j1,j4),s(j5,j6)))*za(j2,j5)**2)
     ./
     .(za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))-
     .(Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))*za(j2,j5)*
     .(za(j1,j5)*za(j2,j3)-za(j1,j2)*za(j3,j5)))/
     .(za(j1,j3)*za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))+
     .(Lsm1_2me(t(j1,j2,j3),t(j2,j3,j4),s(j2,j3),s(j5,j6))*za(j2,j5)*
     .(-(za(j1,j5)*za(j2,j4))+za(j1,j2)*za(j4,j5)))/
     .(za(j1,j4)**2*za(j2,j3)*za(j3,j4)*za(j5,j6))+
     .(Lsm1_2me(t(j1,j2,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*za(j2,j5)*
     .(za(j2,j4)*za(j3,j5)+za(j2,j3)*za(j4,j5)))/
     .(za(j1,j4)*za(j2,j3)*za(j3,j4)**2*za(j5,j6))-
     .(2._dp*L0(-t(j2,j3,j4),-s(j5,j6))*za(j1,j5)*za(j2,j5)*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4)))/
     .(s(j5,j6)*za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))-
     .(2._dp*L0(-t(j2,j3,j4),-s(j2,j3))*za(j2,j5)*za(j4,j5)*zb(j3,j4))/
     .(s(j2,j3)*za(j1,j4)*za(j3,j4)*za(j5,j6))

      return
      end
      
