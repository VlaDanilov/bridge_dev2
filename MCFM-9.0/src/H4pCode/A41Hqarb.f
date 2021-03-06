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
 
c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      function A41Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41Hqarbmppm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A41phiqarbmppm

      A41Hqarbmppm=A41phiqarbmppm(j1,j2,j3,j4,za,zb)
     &            +A41phiqarbmppm(j2,j1,j4,j3,zb,za)
      
      return
      end

      function A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41Hqarbmpmp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A41phiqarbmpmp

      A41Hqarbmpmp=A41phiqarbmpmp(j1,j2,j3,j4,za,zb)
     &            +A41phiqarbmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
