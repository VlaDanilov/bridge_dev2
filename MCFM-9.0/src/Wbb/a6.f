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
 
      function a6(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: a6
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     st is a string to choose between pp,pm or sl
*     implementation of Eq. (3.3) of BDKW hep-ph/9610370
*     character string st can take the value pp,pm or sl
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*2 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: atree,vv,switchyard

      a6=atree(st,j1,j2,j3,j4,j5,j6,za,zb)*vv(st,j1,j2,j3,j4,j5,j6)
      a6=a6+switchyard(st,j1,j2,j3,j4,j5,j6,za,zb)
      
      end

      function switchyard(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: switchyard
c-----switchyard function to direct to pm,pp or st
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*2 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: fpp,fpm,fsl
      if     (st == 'pp') then
           switchyard=fpp(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif (st == 'pm') then
           switchyard=fpm(j1,j2,j3,j4,j5,j6,za,zb) 
      elseif (st == 'sl') then
           switchyard=fsl(j1,j2,j3,j4,j5,j6,za,zb) 
      endif
      end



