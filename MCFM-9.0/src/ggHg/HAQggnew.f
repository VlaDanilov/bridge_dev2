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
 
      subroutine HAQggnew(p1,p2,p3,p4,ampsq,ampsq_ba,ampsq_ab,ampsq_sym)
      implicit none
      include 'types.f'
c--- Note that the ab/ba ordering of the color-ordered matrix elements
c--- returned by this routine is OPOPSITE to that in hqqggdfm.f;
c--- the correct ordering is important for real subtraction terms
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: p1,p2,p3,p4,j1,j2,j3
      complex(dp):: ab(2,2,2),ba(2,2,2)
      real(dp):: ampsq,ampsq_ab,ampsq_ba,ampsq_sym

      call Amplo_AQgg(p1,p2,p3,p4,ab,ba)

c--- calculate the matrix element as the sum of two colour orderings,
c--- plus a colour-suppressed QE.e-_dplike piece which is symmetric
c--- in the ordering of the two gluons
      ampsq_ab=0._dp
      ampsq_ba=0._dp
      ampsq_sym=0._dp
      do j1=1,2
      do j2=1,2
      do j3=1,2
      ampsq_ab=ampsq_ab
     &  +cf*xn**2/2._dp*abs(ab(j1,j2,j3))**2
      ampsq_ba=ampsq_ba
     &  +cf*xn**2/2._dp*abs(ba(j1,j2,j3))**2
      ampsq_sym=ampsq_sym
     &  -cf/2._dp*abs(ab(j1,j2,j3)+ba(j1,j2,j3))**2
      enddo
      enddo
      enddo

      ampsq=ampsq_ab+ampsq_ba+ampsq_sym

      return
      end



