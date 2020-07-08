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
 
c--- Variable declarations and common blocks for producing
c--- histograms with irregular bin widths
      logical:: irregbin(maxhisto)
      integer:: nirreg,irregbin_ptr(maxhisto)
c--- Maximum 10 irregular histograms, with maximum 30 bins in each
      real(dp):: irregbinedges(10,30)
      common/irregbinmcfm/irregbinedges,irregbin,nirreg,irregbin_ptr
      
