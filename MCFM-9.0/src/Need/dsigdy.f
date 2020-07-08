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
 
      function dsigdy(x)      
      implicit none
      include 'types.f'
      include 'nyy.f'
      real(dp)::dsigdy
      real(dp)::sig(nyy),xyy(nyy),err
      logical:: first
      data first/.true./
      save sig,xyy
      if (first) then
      first=.false.
      open(unit=47,file='outw+.dat',status='old')
      do ny=1,nyy
      read(47,*) xyy(ny),sig(ny),err
c      write(6,*) xyy(ny),sig(ny)
      enddo
      close(unit=47)
      endif
      mpot=3
      dsigdy=1d3*ddvdif(sig,xyy,nyy,x,mpot)    
      return       
      end
