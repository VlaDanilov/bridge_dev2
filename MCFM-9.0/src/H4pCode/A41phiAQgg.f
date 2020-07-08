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
 
      function A41phiAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiAQggmpmm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiAQggmpmmL,A1phiAQggmpmmR,A1phiAQggmpmmF

c--- implementation of arXiv:0906.0008, Eq.(2.21)
      A41phiAQggmpmm=
     &     +A1phiAQggmpmmL(j1,j2,j3,j4,za,zb)
     &     -A1phiAQggmpmmR(j1,j2,j3,j4,za,zb)/xn**2
     &     +A1phiAQggmpmmF(j1,j2,j3,j4,za,zb)*real(nflav,dp)/xn

      return
      end
      
      function A41phiAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiAQggmppp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiAQggmpppL,A1phiAQggmpppR,A1phiAQggmpppF

c--- implementation of arXiv:0906.0008, Eq.(2.21)

      A41phiAQggmppp=
     &     +A1phiAQggmpppL(j1,j2,j3,j4,za,zb)
     &     -A1phiAQggmpppR(j1,j2,j3,j4,za,zb)/xn**2
     &     +A1phiAQggmpppF(j1,j2,j3,j4,za,zb)*real(nflav,dp)/xn

      return
      end
      
      function A41phiAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiAQggmpmp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiAQggmpmpL,A1phiAQggmpmpR,A1phiAQggmpmpF

c--- implementation of arXiv:0906.0008, Eq.(2.21)
      A41phiAQggmpmp=
     &     +A1phiAQggmpmpL(j1,j2,j3,j4,za,zb)
     &     -A1phiAQggmpmpR(j1,j2,j3,j4,za,zb)/xn**2
     &     +A1phiAQggmpmpF(j1,j2,j3,j4,za,zb)*real(nflav,dp)/xn

      return
      end
      
      function A41phiAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiAQggmppm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1phiAQggmppmL,A1phiAQggmppmR,A1phiAQggmppmF

c--- implementation of arXiv:0906.0008, Eq.(2.21)
      A41phiAQggmppm=
     &     +A1phiAQggmppmL(j1,j2,j3,j4,za,zb)
     &     -A1phiAQggmppmR(j1,j2,j3,j4,za,zb)/xn**2
     &     +A1phiAQggmppmF(j1,j2,j3,j4,za,zb)*real(nflav,dp)/xn

      return
      end
      
