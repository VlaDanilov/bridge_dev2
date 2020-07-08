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
 
      subroutine amplonumer(amp)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp):: amp(3,2,2,2,2)
      amp(1,2,2,2,2)=
     & + 12618771520986347._dp/21822156628165._dp
     & - 11092695284745761._dp/21822156628165._dp*im
      amp(1,1,1,2,2)=
     & - 3264288934742897417559._dp/17715248572900975165._dp
     & - 435088026732024661529._dp/3543049714580195033._dp*im
      amp(1,1,2,1,2)=
     & - 1370264214._dp/26881165._dp
     & - 97312640._dp/5376233._dp*im
      amp(1,1,2,2,2)=
     & - 1892123283187172788443090._dp/16364303141407230480677._dp
     & + 372907569759498865644745607._dp/654572125656289219227080._dp*im
      amp(1,1,1,1,2)=
     & + 557316427894353929295783568._dp/49726517903897927585527335._dp
     & - 9418892072194906502811319261._dp/745897768558468913782910025._dp*im
      return
      end
