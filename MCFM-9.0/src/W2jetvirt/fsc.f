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
 
      function Fsc(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: Fsc
      integer:: j1,j2,j3,j4,j5,j6
      character*9 st
      complex(dp):: Fsc1,Fsc2,Fsc3,Fsc4,Fsc5,Fsc6,Fsc7,Fsc8
      if(st=='q+g-g+qb-') then
      Fsc=Fsc1(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+g+g-qb-') then
      Fsc=Fsc2(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+g+g+qb-') then
      Fsc=Fsc3(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+g+qb-g-') then
      Fsc=Fsc4(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+g+qb-g+') then
      Fsc=Fsc5(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+qb-g-g+') then
      Fsc=Fsc6(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+qb-g+g-') then
      Fsc=Fsc7(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif(st=='q+qb-g+g+') then
      Fsc=Fsc8(j1,j2,j3,j4,j5,j6,za,zb) 
      endif
      return
      end
