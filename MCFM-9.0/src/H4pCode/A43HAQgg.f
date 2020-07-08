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
 
      function A43HAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmppp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmm
      integer:: j1,j2,j3,j4
      
      A43HAQggmppp= 
c     &   A43phiAQggmppp(j1,j2,j3,j4,za,zb)  ! This term is zero
     &  -A43phiAQggmpmm(j2,j1,j4,j3,zb,za)
      
      return
      end
      
      function A43HAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmpmm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmm
      integer:: j1,j2,j3,j4
      
      A43HAQggmpmm= 
     &    A43phiAQggmpmm(j1,j2,j3,j4,za,zb)
c     &   -A43phiAQggmppp(j2,j1,j4,j3,zb,za) ! This term is zero
      
      return
      end
      
      function A43HAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmpmp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmp
      integer:: j1,j2,j3,j4
      
      A43HAQggmpmp= 
     &    A43phiAQggmpmp(j1,j2,j3,j4,za,zb)
     &   -A43phiAQggmpmp(j2,j1,j4,j3,zb,za)
      
      return
      end
      
      function A43HAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43HAQggmppm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmppm
      integer:: j1,j2,j3,j4
      
      A43HAQggmppm= 
     &    A43phiAQggmppm(j1,j2,j3,j4,za,zb)
     &   -A43phiAQggmppm(j2,j1,j4,j3,zb,za)
      
      return
      end
      
