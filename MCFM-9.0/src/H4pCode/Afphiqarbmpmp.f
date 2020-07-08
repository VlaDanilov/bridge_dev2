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
 
      function Afphiqarbmpmp(j1,j2,j3,j4,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Afphiqarbmpmp
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
C---arXIv:09060008v1 Eq.(4.15)
      integer:: j1,j2,j3,j4
      complex(dp):: Afphiqarbmppm
      Afphiqarbmpmp=-Afphiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

