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
 
      function atreepp(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: atreepp
      
c---Atreepp is the amplitude for
c---q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
c---All outgoing particles are right-handed, except q(p3)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c-----Note all tree amplitudes divided by im wrt to BDK paper
      atreepp=
     & +zb(j1,j2)*za(j5,j4)*(za(j3,j1)*zb(j1,j6)+za(j3,j2)*zb(j2,j6))
     &      /(s(j2,j3)*s(j5,j6)*t(j1,j2,j3))
     & +za(j3,j4)*zb(j6,j1)*(za(j5,j3)*zb(j3,j2)+za(j5,j4)*zb(j4,j2))
     &       /(s(j2,j3)*s(j5,j6)*t(j2,j3,j4))
      return
      end
