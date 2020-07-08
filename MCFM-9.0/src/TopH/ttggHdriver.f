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
 
      subroutine ttggHdriver(q,ampsq1,ampsq2,ampsq0)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: q(mxpart,4),ampsq1,ampsq2,ampsq0
      integer:: h1,h2,h3,h4
      complex(dp)::
     & ampAB(2,2,2,2),ampBA(2,2,2,2),ampQED(2,2,2,2)
      call ttggHamp(q,3,4,1,2,1,1,ampAB)
      call ttggHamp(q,3,4,2,1,1,1,ampBA)
      ampsq1=0d0
      ampsq2=0d0
      ampsq0=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      ampQED(h1,h2,h3,h4)=ampAB(h1,h2,h3,h4)+ampBA(h1,h2,h4,h3)
      ampsq1=ampsq1
     &  +real(ampAB(h1,h2,h3,h4)*conjg(ampAB(h1,h2,h3,h4)))
      ampsq2=ampsq2
     &  +real(ampBA(h1,h2,h4,h3)*conjg(ampBA(h1,h2,h4,h3)))
      ampsq0=ampsq0
     & -real(ampQED(h1,h2,h3,h4)*conjg(ampQED(h1,h2,h3,h4)))/xn**2
      enddo
      enddo
      enddo
      enddo
      return
      end
