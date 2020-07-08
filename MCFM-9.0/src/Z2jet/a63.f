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
 
      function a63(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
C     implementation of last parts of 
C     Eqs. (2.15) and (2.16 ) of Nucl Phys. 513, 3 (1998)
      include 'types.f'
      complex(dp):: a63
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: a6ax
      integer:: j1,j2,j3,j4,j5,j6
      character*2 st

      if (st == 'pp') then
      a63=-a6ax(j1,j4,j3,j2,j5,j6,za,zb)
      elseif (st == 'pm') then
      a63=+a6ax(j1,j4,j2,j3,j5,j6,za,zb)
      else
      write(6,*) 'Unimplemented st in a63',st 
      stop
      endif
      return
      end

