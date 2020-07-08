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
 
      subroutine idjet(pjet,jetindex,numljets,numbjets,
     &     pljet,pbjet)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      real(dp):: pjet(mxpart,4),pljet(mxpart,4),pbjet(mxpart,4)
      integer:: numbjets,numljets,j,jetindex(mxpart)

      numljets=0
      numbjets=0

      do j=1,jets
         if (jetlabel(j) == 'pp' .or. jetlabel(j) == 'qj') then 
            numljets=numljets+1
            pljet(numljets,:)=pjet(jetindex(j),:)
         elseif (jetlabel(j) == 'bq' .or. jetlabel(j) == 'ba'
     &      .or. jetlabel(j) == 'bb') then
            numbjets=numbjets+1
            pbjet(numbjets,:)=pjet(jetindex(j),:)
         else
            write(*,*) 'In idjet, something wrong in jetlabel', 
     &           jetlabel(j)
            stop
         endif
      enddo
      
      return
      end
