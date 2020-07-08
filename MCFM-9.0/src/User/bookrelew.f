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
 
      subroutine bookrelew(n,tag,title,var,wt,wtnoew,xmin,xmax,xbin)
c--- routine to set up histograms that contain the relative effect
c--- of electroweak corrections; should be called from nplotter
      implicit none
      include 'types.f'
      integer n
      integer tag
      character*(*) title
      real(dp):: var,wt,wtnoew,xmin,xmax,xbin,wtwithew
      
      wtwithew=wtnoew+wt
      
      call bookplot(n,tag,trim(title)//' - no EW',var,wtnoew,wtnoew**2,
     &               xmin,xmax,xbin,'lin')
      n=n+1
      call bookplot(n,tag,trim(title)//' - with EW',var,wtwithew,wtwithew**2,
     &               xmin,xmax,xbin,'lin')
      n=n+1
      call bookplot(n,tag,trim(title)//' - +RELEW+',var,wtwithew,wtwithew**2,
     &               xmin,xmax,xbin,'lin')
      n=n+1
      
      return
      end
      
