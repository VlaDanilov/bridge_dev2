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
 
      function Hggggvsq(j1,j2,j3,j4)
      implicit none
      include 'types.f'
      real(dp):: Hggggvsq
      
      include 'epinv.f'
      include 'epinv2.f'
      !include 'first_time.f' 
C---  matrix element squared for H--g(j1)+g(j2)+g(j3)+g(j4)
      integer:: al,j1,j2,j3,j4
      real(dp):: q(4,5),qswap(4,5),sqres(-2:0)
      common/GZmom/q
!$omp threadprivate(/GZmom/)
      
      do al=1,4
      qswap(al,1)=q(al,j1)
      qswap(al,2)=q(al,j2)
      qswap(al,3)=q(al,j3)
      qswap(al,4)=q(al,j4)
      qswap(al,5)=-q(al,j1)-q(al,j2)-q(al,j3)-q(al,j4)
      enddo 

C      write(*,*) 'scale, hmass',scale, hmass
      call GZHggggvsqPoles(qswap,sqres) 
C      call GZHggggvsq(qswap,sqres,first_time,new_event,scale,hmass) 
      Hggggvsq=epinv*epinv2*sqres(-2)+epinv*sqres(-1)+sqres(0)
      !if (first_time) first_time = .false. 
      !if (new_event) new_event = .false. 
      return
      end


