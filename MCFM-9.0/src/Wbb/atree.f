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
 
      function atree(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: atree
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---- the basic process is atreepp which is the amplitude for
c---- q-(p4)+Q(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
c---- All outgoing particles are right-handed, 
c---- except Q(p3) which is left-handed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      character*2 st
      real(dp):: t
      complex(dp):: tree
c-----Note all tree amplitudes divided by im wrt to BDKW paper
c-----Statement function for pp (others derived by crossing)
c-----This is an implementation of BDKW Eq.(3.10)
      tree(j1,j2,j3,j4,j5,j6)=
     & +zb(j1,j2)*za(j5,j4)*(za(j3,j1)*zb(j1,j6)+za(j3,j2)*zb(j2,j6))
     &      /(s(j2,j3)*s(j5,j6)*t(j1,j2,j3))
     & +za(j3,j4)*zb(j6,j1)*(za(j5,j3)*zb(j3,j2)+za(j5,j4)*zb(j4,j2))
     &       /(s(j2,j3)*s(j5,j6)*t(j2,j3,j4))

      if (st == 'pp') then
         atree=+tree(j1,j2,j3,j4,j5,j6)
      elseif (st == 'pm') then
         atree=-tree(j1,j3,j2,j4,j5,j6)
      elseif (st == 'sl') then
c---this is not a real tree amplitude but it is an auxiliary quantity.
         atree=-tree(j3,j1,j2,j4,j5,j6)
      endif

      return 
      end
