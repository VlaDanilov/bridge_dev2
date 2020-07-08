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
 
C      function Hqarbvsqnum(j1,j2,j3,j4)
C      implicit none
      include 'types.f'
      real(dp):: Hqarbvsqnum
C      
C      include 'epinv.f'
C      include 'first_time.f'
CC---  matrix element squared for H--q(j1)+a(j2)+g(j3)+g(j4)
C      integer:: al,j1,j2,j3,j4
C      real(dp):: q(4,5),qswap(4,5),sqres(-2:0)
C      common/GZmom/q
C      do al=1,4
C      qswap(al,1)=q(al,j1)
C      qswap(al,2)=q(al,j2)
C      qswap(al,3)=q(al,j3)
C      qswap(al,4)=q(al,j4)
C      qswap(al,5)=-q(al,j1)-q(al,j2)-q(al,j3)-q(al,j4)
C      enddo
C      call GZHqarbvsq(qswap,sqres,first_time,new_event) 
C      Hqarbvsqnum=epinv**2*sqres(-2)+epinv*sqres(-1)+sqres(0)
C      if (first_time) first_time = .false. 
C      if (new_event) new_event = .false. 
C      return
C      end
