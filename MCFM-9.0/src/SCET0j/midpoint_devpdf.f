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
 
! Simple derivative by Taylor expansion about midpoint:
!    f'(x)= 1/(2a)*(f(a+x)-f(x-a)) + O(f''')
      subroutine midpoint_devpdf(ih1,ibeam,x,pdf)
      implicit none
      include 'types.f'
      include 'scale.f'
      include 'facscale.f'
      real(dp) :: x,pdf(-5:5),upx,downx
      real(dp) :: uppdfs(-5:5),downpdfs(-5:5)
      real(dp) :: smallxfac,tinyx
      integer ih1,ibeam,j

      smallxfac=1e-4_dp

      tinyx=smallxfac*x
      upx=x+tinyx
      downx=x-tinyx
      
      call fdist(ih1,upx,facscale,uppdfs)
      call fdist(ih1,downx,facscale,downpdfs)

      pdf(:)=0._dp
      do j=-5,5
         pdf(j)=(uppdfs(j)-downpdfs(j))/(2._dp*tinyx)
      enddo

      return
      end
