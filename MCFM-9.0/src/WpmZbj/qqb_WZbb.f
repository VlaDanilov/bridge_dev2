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
 
      subroutine qqb_WZbb(p,msq)
      implicit none
      include 'types.f'
      
!---  Author: R.K. Ellis February 2013
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'nwz.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),WZbbmsq

      call spinoru(8,p,za,zb)

c--- initialize matrix elements to zero
      msq(:,:)=0._dp

      if (nwz == 1) then
      msq(+2,-1)=WZbbmsq(1,2,3,4,5,6,7,8)
      msq(-1,+2)=WZbbmsq(2,1,3,4,5,6,7,8)

      msq(+4,-3)=msq(+2,-1)
      msq(-3,+4)=msq(-1,+2)

      elseif (nwz == -1) then

      msq(+1,-2)=WZbbmsq(1,2,3,4,5,6,7,8)
      msq(-2,+1)=WZbbmsq(2,1,3,4,5,6,7,8)

      msq(+3,-4)=msq(+1,-2)
      msq(-4,+3)=msq(-2,+1)

      endif

      return
      end


