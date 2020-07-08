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
 
      subroutine storecsz(mcs)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the Z2jet matrix elements into separate arrays for each
c-- incoming parton case
      
      include 'mmsq_cs.f'
      integer:: i,j,k
      real(dp):: mcs(0:2,2,2)
      
      do i=0,2
        do j=1,2
          do k=1,2
        mcs(i,j,k)=mmsq_cs(i,j,k)
          enddo
        enddo
      enddo
      return
      end
