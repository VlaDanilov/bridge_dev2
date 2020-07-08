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
 
c-----multiplication of a barred spinor with k-slash from the right
C-----and return resultant spinor f. Weyl representation.
C     Energy component in MCFM notation = k(4)
      subroutine Ubkslash(spinor,k,f) 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'swapxz.f'
      include 'debug.f'
      complex(dp):: spinor(4),k(4),f(4),kslash(4,4),E,kx,ky,kz
      integer:: i,j
      logical,save::first
      data first/.true./
!$omp threadprivate(first)

      if (first) then
        if (debug) write(6,*) 'Ubkslash:swapxz=',swapxz
        first=.false.
      endif

      if (swapxz) then
C----create kslash after performing the swap (x<->z),(y->-y)
      E=k(4)
      kx=+k(3)
      ky=-k(2)
      kz=+k(1)
      else
      E=k(4)
      kx=+k(1)
      ky=+k(2)
      kz=+k(3)
      endif

      kslash(1,1)=czip
      kslash(1,2)=czip
      kslash(1,3)=E+kz
      kslash(1,4)=kx-im*ky

      kslash(2,1)=czip
      kslash(2,2)=czip
      kslash(2,3)=kx+im*ky
      kslash(2,4)=E-kz

      kslash(3,1)=E-kz
      kslash(3,2)=-kx+im*ky
      kslash(3,3)=czip
      kslash(3,4)=czip

      kslash(4,1)=-kx-im*ky
      kslash(4,2)=E+kz
      kslash(4,3)=czip
      kslash(4,4)=czip

      do i=1,4
      f(i)=czip
      do j=1,4
      f(i)=f(i)+spinor(j)*kslash(j,i)
      enddo
      enddo  
      return
      end



