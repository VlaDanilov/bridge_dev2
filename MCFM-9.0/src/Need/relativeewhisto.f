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
 
      subroutine relativeewhisto(N)
c--- Compute relative EW corrections in histogram N,
c--- assuming (NOEW+EW) is currently in N and (NOEW) is in N-2
      implicit none
      include 'types.f'
      include 'histo.f'
      integer N,j
      double precision xrel,xrelsig
      
      do j=1,NBIN(N)
      xrel=HIST(N,j)/HIST(N-2,j)-1d0
c--- at this point, histogram errors are stored in 2*maxhisto+N
      xrelsig=
     & xrel*sqrt((HIST(2*maxhisto+N,j)/HIST(N,j))**2
     &          +(HIST(2*maxhisto+N-2,j)/HIST(N-2,j))**2)
      HIST(N,j)=xrel
      HIST(2*maxhisto+N,j)=abs(xrelsig)
      enddo
      
      j=index(title(N),'+RELEW+')
      title(N)(j:j+6)='rel EW '
      
      return
      end
      
      
