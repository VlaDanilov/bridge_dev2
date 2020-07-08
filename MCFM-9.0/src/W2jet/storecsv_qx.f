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
 
      subroutine storecsv_qx(i,j)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the W2jet_gvec matrix elements into elements (..,i,j) of q1q2
      
      include 'mmsqv_cs.f'
      integer:: i,j,k
      real(dp):: q1q2(0:2,-1:1,-1:1)
      common/q1q2/q1q2
!$omp threadprivate(/q1q2/)
      
      do k=0,2
        q1q2(k,i,j)=mmsqv_cs(k,+1,+1)
      enddo
      
      return
      end
