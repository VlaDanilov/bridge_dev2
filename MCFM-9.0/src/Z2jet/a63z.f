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
 
      function a63z(i1,i2,i3,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a63z
c--- this function should be passed the (true) helicities of the particles,
c--- namely with i1 incoming and i2,i3 outgoing (1=left, 2=right)
c--- to compare with BDKW paper, first flip i1 to outgoing and then
c--- note that left corresponds to + and right to -
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,nhel
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a63

      nhel=(i1-1)+(i2-1)*2+(i3-1)*4

c--- note that reversing all helicities swaps za,zb
c--- and switching 3rd helicity swaps j5 and j6      
      if     (nhel == 0) then
        a63z=+a63('pm',j1,j2,j3,j4,j5,j6,zb,za)
      elseif (nhel == 1) then
        a63z=+a63('pp',j1,j2,j3,j4,j6,j5,za,zb)
      elseif (nhel == 2) then
        a63z=+a63('pp',j1,j2,j3,j4,j5,j6,zb,za)
      elseif (nhel == 3) then
        a63z=+a63('pm',j1,j2,j3,j4,j6,j5,za,zb)
      elseif (nhel == 4) then
        a63z=+a63('pm',j1,j2,j3,j4,j6,j5,zb,za)      
      elseif (nhel == 5) then
        a63z=+a63('pp',j1,j2,j3,j4,j5,j6,za,zb)
      elseif (nhel == 6) then
        a63z=+a63('pp',j1,j2,j3,j4,j6,j5,zb,za)
      elseif (nhel == 7) then
        a63z=+a63('pm',j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(*,*) 'Unsupported helicity assignment'
        stop
      endif

      return
      end
     
