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
 
      function ampqqbgll(p1,p2,p3,p4,p5,al,be,ga,de,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
      integer:: p1,p2,p3,p4,p5
      complex(dp)::ampqqbgll
      complex(dp)::al,be,ga,de
      ampqqbgll=
     & +al*za(p1,p2)*zb(p1,p5)*za(p4,p2)/(za(p1,p3)*za(p3,p2))
     & +be*zb(p3,p5)*za(p4,p2)/za(p1,p3)
     & +ga*zb(p3,p1)*zb(p3,p5)*za(p4,p1)/(za(p1,p3)*zb(p2,p3))
     & +de*zb(p3,p1)/(za(p1,p3)*zb(p2,p1))
     & *(zb(p1,p5)*za(p4,p1)+zb(p2,p5)*za(p4,p2)+zb(p3,p5)*za(p4,p3))
     
      return
      end
