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
 
      subroutine arraysort(imax,vals,sortorder)
      implicit none
      include 'types.f'
c--- Given a a 1-dimensional array of quantities, vals(1:mxpart),
c--- this routine sorts the first imax entries such that
c--- vals(sortorder(1)) > vals(sortorder(2)) >...> vals(sortorder(imax)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: imax,sortorder(mxpart),i1,i2,iswap
      real(dp):: vals(mxpart)

c-- initial ordering      
      do i1=1,imax
      sortorder(i1)=i1
      enddo

c--- catch trivial case      
      if (imax == 1) return
      
      do i1=imax-1,1,-1
      do i2=1,i1
        if (vals(sortorder(i2)) < vals(sortorder(i2+1))) then
          iswap=sortorder(i2)
          sortorder(i2)=sortorder(i2+1)
          sortorder(i2+1)=iswap
        endif
      enddo
      enddo
      
      return
      end
      
