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
 
      function A0phiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiqarbmppm
      
C----Expresssion of Eq. (A.4) of arXiv:0906008
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phiqarbmppm=-za(j1,j4)**2/(za(j1,j2)*za(j3,j4))
      return
      end


      function A0phiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiqarbmpmp
      
C----Expresssion of Eq. (A.4) of arXiv:0906008
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiqarbmppm
      A0phiqarbmpmp=-A0phiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

      function A0phidqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidqarbmppm
      
C----Expresssion of Eq. (A.4) of arXiv:0906008
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phidqarbmppm=-zb(j2,j3)**2/(zb(j1,j2)*zb(j3,j4))
      return
      end

      function A0phidqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidqarbmpmp
      
C----Expresssion of Eq. (A.4) of arXiv:0906008
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A0phidqarbmppm
      integer:: j1,j2,j3,j4
      A0phidqarbmpmp=-A0phidqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end
