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
 
      function A41phiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiqarbmppm
C     arXiv:09060008v1, Eq.2.14
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'zprods_decl.f'
      complex(dp):: Alcphiqarbmppm,Alcphiqarbmpmp,Aslcphiqarbmppm,
     & Afphiqarbmppm
      integer:: j1,j2,j3,j4
      A41phiqarbmppm=
     &                      Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &            -2._dp/xnsq*(Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &                      +Alcphiqarbmpmp(j1,j2,j4,j3,za,zb))
     &            -1._dp/xnsq*Aslcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &    +real(nflav,dp)/xn*Afphiqarbmppm(j1,j2,j3,j4,za,zb)
      return
      end

      function A41phiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiqarbmpmp
C     arXiv:09060008v1, Eq.2.15
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'zprods_decl.f'
      complex(dp):: Alcphiqarbmpmp,Alcphiqarbmppm,Aslcphiqarbmppm,
     & Afphiqarbmpmp
      integer:: j1,j2,j3,j4
      A41phiqarbmpmp=
     &                      Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &            -2._dp/xnsq*(Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &                      +Alcphiqarbmppm(j1,j2,j4,j3,za,zb))
     &            +1._dp/xnsq*Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
     &    +real(nflav,dp)/xn*Afphiqarbmpmp(j1,j2,j3,j4,za,zb)
      return
      end
