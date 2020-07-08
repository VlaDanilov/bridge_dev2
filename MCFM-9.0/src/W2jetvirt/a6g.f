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
 
      function a6g(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: a6g
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     August, 1999.                                                    *
*     implementation of Eq. (6.1) of BDKW hep-ph/9708239               *
*     character string st can take various values                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*9 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeg,vvg,fcc,fsc

      a6g=a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)*vvg(st,j1,j2,j3,j4,j5,j6)
      a6g=a6g+fcc(st,j1,j2,j3,j4,j5,j6,za,zb)
     &       +fsc(st,j1,j2,j3,j4,j5,j6,za,zb)

      end

