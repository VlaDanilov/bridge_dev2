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
 
      function F31m(s)
      implicit none
      include 'types.f'
      complex(dp):: F31m
      
      include 'epinv.f'
      include 'scale.f'
      real(dp):: s
c      complex(dp):: qlI3
      complex(dp):: lnrat
c      integer:: ep
c      F31m=czip
c      do ep=-2,0
c      F31m=F31m+s*epinv**(-ep)*qlI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,ep)
c      enddo
      
c--- NOTE: checked on 8/31/09 that this agrees with the expression above
      F31m=epinv**2-epinv*lnrat(-s,musq)+0.5_dp*lnrat(-s,musq)**2
      
      return
      end

