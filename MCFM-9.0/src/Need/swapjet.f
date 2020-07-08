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
 
      subroutine swapjet(pjet,jetindex,i,j)
        use jettagging
      implicit none
c--- swaps jets i..j in pjet
      
      include 'constants.f'
      include 'nf.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,j,k,jetindex(mxpart)
      real(dp):: pjet(mxpart,4),tmp
      integer :: itmp
      character*2 chartmp
 
c--- escape if we're trying to swap the same jets
      if (i == j) return

      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
 
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp

      itmp = jetcontent(i)
      jetcontent(i) = jetcontent(j)
      jetcontent(j) = itmp
      
      return
      end
