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
 
      function mcfmccdb0(psq,m1sq,m2sq)
! Implementation of bubble derivative: d/dpsq B0(psq,m1sq,m2sq)
! from Eq. (4.25) of A. Denner, arXiv:0709.1075
      implicit none
      include 'types.f'
      real(dp):: psq,m1sq,m2sq
      complex(dp):: mcfmccdb0,bit,root,rp,rm,r
      real(dp), parameter:: tiny=1.e-10_dp

!--- work out the value of r from the defining equation:
!     r + 1/r = (m1sq+m2sq-psq-ie)/(m0*m1)
      bit=(m1sq+m2sq-psq-cmplx(0._dp,1.e-5_dp,dp))/sqrt(m1sq*m2sq)
      
      root=sqrt(bit**2-4._dp)
      rp=(bit+root)/2._dp
      rm=bit-rp

! For psq=0, implement special case of m1sq=m2sq
      if (abs(psq) < tiny) then
        if (abs(m1sq-m2sq) < tiny) then
          mcfmccdb0=cmplx(1._dp/(6._dp*m1sq),0._dp,dp)
          return
        else
          write(6,*) 'mcfmccdb0 not implemented for psq=0, m1sq /= m2sq'
          stop
        endif
      endif

! choose root -- could be done with more attention to numerical precision ...
      r=rm

      mcfmccdb0=-(m1sq-m2sq)/psq**2*log(m2sq/m1sq)/2._dp
     & +sqrt(m1sq*m2sq)/psq**2*(1._dp/r-r)*log(r)
     & -(1._dp+(r**2+1._dp)/(r**2-1._dp)*log(r))/psq

      end function mcfmccdb0
