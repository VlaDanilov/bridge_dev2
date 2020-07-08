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
 
      function A42phiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A42phiqarbmppm
C     arXiv:09060008v1, Eq.2.16
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'zprods_decl.f'
      complex(dp):: Alcphiqarbmppm,Alcphiqarbmpmp,Aslcphiqarbmppm,
     & Afphiqarbmppm
      integer:: j1,j2,j3,j4
      A42phiqarbmppm=
     &           Alcphiqarbmpmp(j1,j2,j4,j3,za,zb)
     & +1._dp/xnsq*(Alcphiqarbmpmp(j1,j2,j4,j3,za,zb)
     &           +Alcphiqarbmppm(j1,j2,j3,j4,za,zb))
     & +1._dp/xnsq*Aslcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &    -real(nflav,dp)/xn*Afphiqarbmppm(j1,j2,j3,j4,za,zb)
      return
      end

      function A42phiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A42phiqarbmpmp
C     based on arXiv:09060008v1, Eq.2.17
c       but with slight modifications traced back to equation 2.8 of
c        Z. Bern, L. Dixon, D. Kosower & S. Weinzierl, hep-ph/9610370
c       
c        a) interchanged 3- and 4+ in first term
c        b) changed sign of second line [ i.e. +1/N^2*(Alc+Alc)]

       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'zprods_decl.f'
      complex(dp):: Alcphiqarbmpmp,Alcphiqarbmppm,Aslcphiqarbmppm,
     & Afphiqarbmpmp
      integer:: j1,j2,j3,j4
      A42phiqarbmpmp=
     &           Alcphiqarbmppm(j1,j2,j4,j3,za,zb)
     & +1._dp/xnsq*(Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &           +Alcphiqarbmppm(j1,j2,j4,j3,za,zb))
     & -1._dp/xnsq*Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
     &    -real(nflav,dp)/xn*Afphiqarbmpmp(j1,j2,j3,j4,za,zb)
      return
      end

