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
 
      function F33m(p1sq,p2sq,p3sq)
      implicit none
      include 'types.f'
      complex(dp):: F33m
      
c      include 'scale.f'
      real(dp):: p1sq,p2sq,p3sq
c      complex(dp):: qlI3
      complex(dp):: I3m

c--- NOTE: checked on 8/30/09 that qlI3 == -I3m
c---       and F33m is defined to be (-1)*(scalar integral)
c      F33m=-qlI3(p1sq,p2sq,p3sq,0._dp,0._dp,0._dp,musq,0)
      F33m=I3m(p1sq,p2sq,p3sq)

      return
      end

